# 3MF/STL 广告贴自动分盘脚本

这个仓库包含一个用于把超出 P1S（或任意 25 cm×25 cm 级别打印床）的 3MF **或 STL** 文件自动切成多个可打印托盘的小工具。它的工作方式是：

1. 读取原始 3MF/STL 并提取 Z 方向朝上的表面，把整块广告贴视为在 XY 平面上的 2D 图形；
2. 按打印床尺寸（默认 250 mm×250 mm，可自定义）对 2D 轮廓进行网格化分割，支持自定义相邻分块的重叠量；
3. 对每个分块重新居中到打印床、按原始厚度重新拉伸（extrude）成 3D 模型，然后导出新的 3MF 或 STL 文件。

> ⚠️ **使用前提**：模型需要是类似 0.4 mm 厚的平面薄片，且法向指向 Z 轴。如果是有明显高度变化的立体结构，此脚本不适用。

## 快速开始

1. **安装依赖**

   ```powershell
   python -m venv .venv
   .\.venv\Scripts\activate
   python -m pip install -r requirements.txt
   ```

   如果环境没有 Python，请先从 [python.org](https://www.python.org/downloads/) 安装 3.10+ 版本。

2. **运行脚本**

   ```powershell
   python split_3mf.py your_model.3mf `
       --plate-width-mm 250 `
       --plate-height-mm 250 `
       --overlap-mm 3 `
       --bed-margin-mm 5
   ```

   运行后会在 `your_model_plates/` 目录下生成：

   - `*_rXX_cYY.3mf` / `.stl`：每个分块对应的文件（r=行号，c=列号，方便按原图重新拼装）；
   - `manifest.json`：记录原始文件信息、每个分块所在的全局坐标、行列索引等，可在打印后参考拼贴顺序。

## 常用参数

| 参数 | 默认值 | 说明 |
| --- | --- | --- |
| `--plate-width-mm` / `--plate-height-mm` | 250 | 单个打印床可用尺寸，单位 mm。 |
| `--overlap-mm` | 2 | 分块与相邻分块之间的重叠宽度（便于贴合时切齐）。 |
| `--bed-margin-mm` | 3 | 将分块放到托盘中央时预留的边距。 |
| `--min-feature-area-mm2` | 4 | 过滤掉比该面积更小的碎片。 |
| `--min-normal-z` | 0.6 | 提取 2D 轮廓时，只有法向 Z 分量 ≥ 该值的面会被视为“上表面”。 |
| `--output-format` | 与输入一致（找不到则为 3MF） | 输出文件的格式，支持 3MF / STL。比如输入 STL 时希望仍旧输出 STL，可指定 `--output-format stl`。 |
| `--dry-run` | — | 只生成 manifest 而不导出模型，用于检查分块数量。 |

## 工作流程

1. 加载 3MF/STL 并合并所有几何体，读取模型厚度；
2. 通过截取所有朝上的三角面并投影到 XY 平面，获得 2D 轮廓；
3. 依据床尺寸 + 重叠量计算行列数，对轮廓做网格裁剪；
4. 将裁剪结果平移到打印床中央，恢复原始厚度，导出 3MF/STL；
5. 生成一个 `manifest.json` 记录分块地图。

## 限制与建议

- 模型必须基本平整（类似 0.4 mm 薄片）。大量起伏或厚度变化的模型无法正确恢复；
- 无论来源是 3MF 还是 STL，脚本都会重新根据 2D 轮廓挤出几何，因此原始模型中的贴图/配色信息不会保留；
- 原始 3MF 如果带有多个独立部件，脚本会把它们视作一个整体的 2D 投影进行切分；
- 如果遇到 “Tile does not fit” 提示，说明分块加上边距仍然超出床尺寸，请降低 `--overlap-mm` 或 `--bed-margin-mm`；
- 生成的模型只包含几何数据，需要重新在 Bambu Studio 中设置打印参数（层高 0.4 mm、颜色等）。

如需进一步定制（例如导出 STL/3MF 以外的格式、控制分块顺序、生成预览 SVG），可在 `split_3mf.py` 中扩展对应函数。祝打印顺利! 🎉
