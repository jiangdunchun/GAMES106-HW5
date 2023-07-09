# QEM

> 1. 计算所有顶点的quadric $Q(\pmb{A},\pmb{b},c)$
> 2. 选出所有有效pair $(\pmb{v_1}, \pmb{v_2})$，计算合并后的 $Q(\overline{\pmb{v}})$，添加到堆中
> 3. 合并 $Q$ 最小的pair $(\pmb{v_1}, \pmb{v_2})$ 至 $\overline{\pmb{v}}$
> 4. 在堆中移除所有无效的pair，将新增的pair添加到堆中；重复步骤3直至网格数量减少到目标网格数