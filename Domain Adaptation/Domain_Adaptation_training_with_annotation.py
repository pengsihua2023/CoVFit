#!/usr/bin/env python3
# coding: utf-8
# 第一行是shebang行，它告诉系统这个脚本应该用哪个解释器执行，这里是Python 3。第二行指定了文件的编码格式为UTF-8。

from datasets import Dataset
# 引入Hugging Face的datasets库，用于数据的加载和处理。

import pandas as pd
import numpy as np
import torch
# 引入Pandas、NumPy和PyTorch库，这些是进行数据处理和机器学习的常用库。

from torch.optim import AdamW
# 从PyTorch库中引入AdamW优化器，常用于深度学习模型的参数优化。

from transformers import AutoTokenizer, DataCollatorForLanguageModeling, EsmForMaskedLM, get_cosine_schedule_with_warmup, Trainer, TrainingArguments
# 引入transformers库中的相关类和函数，用于文本处理、模型训练和调度。

from Bio import SeqIO as sio
# 从Biopython库引入SeqIO模块，用于读取和处理生物序列数据。

import random
import math
import datetime
# 引入random、math和datetime库，用于生成随机数、数学运算和处理日期时间。

### Parameters
# 定义训练模型相关的参数。

version = 'ESM2_650M_DA_combo_e30_lrsWarmupCosine'
model_type = 'facebook/esm2_t33_650M_UR50D'
max_length = 1000
batch_size = 5
warmup_epochs = 1
num_epochs = 30
randseed = 13
# 模型版本、类型、最大输入长度、批次大小、预热周期数、训练周期数和随机种子。

### Setting random seeds for reproducibility
# 设置随机种子以保证实验可重复性。

torch.backends.cudnn.deterministic = True
# 使PyTorch在使用CUDA时的行为更加确定性。

random.seed(randseed)
torch.manual_seed(randseed)
torch.cuda.manual_seed(randseed)
np.random.seed(randseed)
# 设置Python内置随机库、PyTorch和NumPy的随机种子。

### Data Prep
# 数据准备部分。

df_cv = pd.read_csv('df_coronaviradae_new.csv')  # coronaviradae sequences
# 读取冠状病毒科的序列数据。

df_sc2 = pd.DataFrame.from_dict({i.id:str(i.seq) for i in sio.parse('legacy_seq_20220831_db99.fasta', 'fasta')}, orient='index', columns=['seq'])  # SARS-CoV-2 sequences
# 从FASTA文件中读取SARS-CoV-2的序列，并转换为DataFrame。

df_uniref = pd.read_table('../uniref/uniref_subset_seq.tsv')  # uniref50 sequences
# 读取UniRef50数据库的序列数据。

# Combine cornoviradae and SARS-CoV-2 sequences to make training set
X_train = df_cv.seq.tolist() + df_sc2.seq.tolist()
train_df = pd.DataFrame(X_train)
train_df.columns = ['seq']
# 将冠状病毒科和SARS-CoV-2的序列合并成训练集。

# Load main dataset and split data by date to make an evaluation dataset of "future" sequences
full_df = pd.read_table('metadata.representative.all_countries.with_date.with_aligned_seq.with_immune_escape_score.txt')
input_df = full_df[['relative_Re', 'country', 'seq', 'date.first', 'relative_Re_log', 'clade']]
cutoff = pd.to_datetime('2022-09-01')
input_df['date.first'] = pd.to_datetime(input_df['date.first'])
recent_df = input_df[input_df['date.first'] >= cutoff]
# 加载包含时间信息的主数据集，并根据日期切分数据，以构建一个包含“未来”序列的评估数据集。

# Input data into dataset.Dataset to prevent incompatibility with slow tokenizer and data_collator
train_dataset = Dataset.from_pandas(train_df[['seq']])
recent_dataset = Dataset.from_pandas(recent_df[['seq']])
test_dataset = Dataset.from_pandas(df_uniref[['seq']])
# 将数据转换为Hugging Face的Dataset格式，以避免与慢速分词器和数据整理工具的不兼容问题。

# Tokenize the datasets
tokenizer = AutoTokenizer.from_pretrained(model_type)
# 使用预训练的分词器初始化。

def tokenizer_for_map(input):  # Tokenizer and params including special_tokens_mask required for MLM
    return tokenizer(
        input['seq'],
        padding='max_length',
        truncation=True,
        max_length=1000,
        return_special_tokens_mask=True
    )
# 定义一个函数用于映射分词过程，包括填充、截断和返回特殊令牌掩码。

column_names = train_dataset.column_names  # This will be the names of all the old columns, to then be deleted after the new tokenized columns are added.
# 获取数据集中的列名，这些列在添加了新的分词列后将被删除。

train_dataset = train_dataset.map(
    tokenizer_for_map,
    batched=True,
    num_proc=8,
    remove_columns=column_names,
)
recent_dataset = recent_dataset.map(
    tokenizer_for_map,
    batched=True,
    num_proc=8,
    remove_columns=column_names,
)
test_dataset = test_dataset.map(
    tokenizer_for_map,
    batched=True,
    num_proc=8,
    remove_columns=column_names,
)
# 使用分词函数映射到训练、最近和测试数据集上，并移除原始列。

### Training
# 训练部分。

data_collator = DataCollatorForLanguageModeling(tokenizer=tokenizer, return_tensors='pt', mlm_probability=0.15)  # Provides random masking and returns tensors during training per-batch
# 初始化数据整理器，用于在训练过程中对输入数据进行随机掩码和批处理张量化。

model = EsmForMaskedLM.from_pretrained(model_type)  # Model with MaskedLM layer on top to provide predictions of masked sequences
# 从预训练模型加载带有掩码语言模型层的模型，用于预测被掩码的序列。

params = filter(lambda x: x.requires_grad, model.parameters())
optimizer = torch.optim.AdamW(params, lr=2e-5)
warmup_steps = math.ceil((train_dataset.num_rows / batch_size) * warmup_epochs)
training_steps = math.ceil((train_dataset.num_rows / batch_size) * num_epochs)
scheduler = get_cosine_schedule_with_warmup(optimizer, num_warmup_steps=warmup_steps, num_training_steps=training_steps)
# 设置优化器和调度器，包括学习率、预热步骤和训练步骤的计算。

training_args = TrainingArguments(
    output_dir=f"{version}/trainer",
    overwrite_output_dir=True,
    num_train_epochs=num_epochs,
    per_device_train_batch_size=batch_size,
    evaluation_strategy='no',
    do_eval=False,
    save_total_limit=1,
    seed=randseed,
    dataloader_num_workers=4,
)
# 定义训练过程的配置参数，包括输出目录、批量大小、评估策略等。

trainer = Trainer(
    model=model,
    args=training_args,
    train_dataset=train_dataset,
    data_collator=data_collator,
    optimizers=[optimizer, scheduler]
)
# 创建训练器实例，并开始训练。

trainer.train()
# 执行训练过程。

trainer.save_model(f"{version}/model")
# 保存训练后的模型。

### Evaluation on uniref50
# 在UniRef50数据集上进行评估。

model = EsmForMaskedLM.from_pretrained(model_type)
# 重新加载原始的预训练模型。

evaluator = Trainer(
  model=model,
  data_collator=data_collator,
  eval_dataset=test_dataset,
  )
# 创建用于评估的训练器实例。

eval_results_old = evaluator.evaluate()
# 执行评估并获取结果。

model = EsmForMaskedLM.from_pretrained(f"{version}/model")
# 加载训练后的模型。

evaluator = Trainer(
  model=model,
  data_collator=data_collator,
  eval_dataset=test_dataset,
  )
# 创建用于评估的训练器实例。

eval_results_new = evaluator.evaluate()
# 执行评估并获取新的结果。

import math
print('ESM2 Uniref Data:')
print(f"Original Model's Perplexity: {math.exp(eval_results_old['eval_loss']):.4f}")
print(f"Adapted Model's Perplexity: {math.exp(eval_results_new['eval_loss']):.4f}")
# 打印UniRef50数据集上的原始模型和训练后模型的困惑度。

### Evaluation on "future" SARS-CoV-2 sequences
# 在“未来”SARS-CoV-2序列数据集上进行评估。

model = EsmForMaskedLM.from_pretrained(model_type)
# 重新加载原始的预训练模型。

evaluator = Trainer(
  model=model,
  data_collator=data_collator,
  eval_dataset=recent_dataset,
  )
# 创建用于评估的训练器实例。

recent_results_old = evaluator.evaluate()
# 执行评估并获取结果。

model = EsmForMaskedLM.from_pretrained(f"{version}/model")
# 加载训练后的模型。

evaluator = Trainer(
  model=model,
  data_collator=data_collator,
  eval_dataset=recent_dataset,
  )
# 创建用于评估的训练器实例。

recent_results_new = evaluator.evaluate()
# 执行评估并获取新的结果。

print('Future SARS-CoV-2 Data:')
print(f"Original Model's Perplexity: {math.exp(recent_results_old['eval_loss']):.4f}")
print(f"Adapted Model's Perplexity: {math.exp(recent_results_new['eval_loss']):.4f}")
# 打印“未来”SARS-CoV-2序列数据集上的原始模型和训练后模型的困惑度。
