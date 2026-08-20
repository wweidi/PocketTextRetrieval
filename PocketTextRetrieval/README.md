# PocketTextRetrieval

口袋级文本–结构双向检索实验项目。

## 当前数据策略

- 主训练集：BioLiP 的实验蛋白–配体口袋，结合配体和局部结合残基构造口袋级文本；不使用蛋白质整体功能描述。
- 主验证/测试集：BioLiP 按 UniProt 分组得到的 pocket-text 验证集和测试集。
- 外部测试集：Receptor.AI 文献口袋数据；OneProt 官方测试集单独用于原始蛋白文本–口袋兼容性比较。
- OneProt 口袋级子集：只有能够将 OneProt 口袋可靠映射到局部结合注释的样本才进入该子集。
- 结构预训练：PLINDER/UniSite-DS 暂作为可选来源，不在第一阶段下载完整数据，避免占用不必要的磁盘空间。

## 目录

```text
configs/       实验配置
data/raw/      原始下载文件
data/interim/  中间映射和清洗文件
data/processed/最终训练/验证/测试索引
data/manifests/数据版本和文件清单
scripts/       下载、映射、清洗和检查脚本
src/           后续模型与数据集代码
logs/          下载和处理日志
results/       实验结果
docs/          数据协议与实验记录
```

## 预训练权重

三套编码器权重下载脚本位于 `scripts/`。下载任务运行在 tmux 会话
`pocket_weight_downloads` 中，支持中断后重新运行并从 `.part` 文件继续下载：

```bash
bash scripts/start_model_weight_downloads.sh
/mnt/HDD0/home/zf25/miniconda3/envs/protein/bin/python scripts/check_model_weights.py
watch -n 10 '/mnt/HDD0/home/zf25/miniconda3/envs/protein/bin/python scripts/check_model_weights.py'
tail -f logs/model_weights_profsa.log
tmux attach -t pocket_weight_downloads
```

目标文件为 ProFSA 的 `profsa_last.ckpt`、Uni-Mol 的
`pocket_pre_220816.pt`，以及 OneProt 使用的 BiomedBERT
`pytorch_model.bin`。

## 数据原则

1. 不把OneProt原始蛋白级文本直接当作口袋级文本。
2. 口袋级文本保留来源和置信度字段。
3. 训练/验证/测试按UniProt蛋白或序列聚类分组，避免同源结构泄漏。
4. OneProt测试蛋白及其高同源蛋白不得进入主训练集。
