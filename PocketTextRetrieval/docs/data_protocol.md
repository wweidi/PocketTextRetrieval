# 数据协议

## 正样本

每条样本必须包含一个文本和一个明确的口袋结构：

```text
text_id, pocket_id, uniprot_id, pdb_id, chain_id,
pocket_residues, pocket_atoms, text, text_source, confidence
```

## 文本版本

- `text_semantic`：不包含精确残基编号，用于主要检索任务。
- `text_full`：可以包含残基、配体和证据细节，用于辅助任务或分析。

## 数据划分

- 训练/验证/测试先按UniProt–PDB二部图的连通分量分组，避免同一蛋白或同一PDB结构跨集合。
- 后续再用MMseqs2生成序列聚类分组，防止高同源蛋白跨集合泄漏。
- OneProt测试蛋白及其高同源蛋白从主训练集排除。
- OneProt-compatible测试保留OneProt原始候选池和匹配关系。
- Gold口袋级测试只保留有可靠口袋文本证据的样本。

## 口袋结构定义

第一阶段使用BioLiP标注的binding-site残基作为口袋掩码，从链级受体结构中提取对应残基；后续可用配体重原子距离阈值生成几何口袋版本。ProFSA和Uni-Mol必须使用完全相同的口袋文件和预处理规则。
