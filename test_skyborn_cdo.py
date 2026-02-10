"""
skyborn-cdo 本地测试脚本
========================
在 skyborn_dev 环境中运行:
    conda activate skyborn_dev
    set PATH=C:\\msys64\\mingw64\\bin;%PATH%
    python test_skyborn_cdo.py
"""

import os
import sys
import glob
import tempfile

# ─── 配置 ───────────────────────────────────────────────
# 确保 MinGW DLL 在 PATH 中（本地开发时需要）
MINGW_BIN = r"C:\msys64\mingw64\bin"
if os.path.isdir(MINGW_BIN) and MINGW_BIN not in os.environ.get("PATH", ""):
    os.environ["PATH"] = MINGW_BIN + os.pathsep + os.environ["PATH"]

GPM_DATA_DIR = r"I:\GPM_Month\Data"
OUTPUT_DIR = r"I:\GPM_Month"


def separator(title: str):
    print(f"\n{'='*60}")
    print(f"  {title}")
    print(f"{'='*60}")


def test_import():
    """测试 1: 包导入"""
    separator("测试 1: 导入 skyborn_cdo")
    from skyborn_cdo import Cdo
    cdo = Cdo()
    print(f"  ✓ 导入成功: {repr(cdo)}")
    return cdo


def test_version(cdo):
    """测试 2: CDO 版本"""
    separator("测试 2: CDO 版本信息")
    version = cdo.version()
    print(f"  版本输出:\n{version}")
    assert "2.5.3" in version, f"版本不匹配: {version}"
    print(f"  ✓ CDO 2.5.3 确认")


def test_operators(cdo):
    """测试 3: 算子列表"""
    separator("测试 3: 算子列表")
    ops = cdo.operators()
    print(f"  可用算子数量: {len(ops)}")
    key_ops = ["sinfo", "mergetime", "remapbil", "sellonlatbox",
               "timavg", "yearmonmean", "fldmean", "griddes"]
    for op in key_ops:
        status = "✓" if op in ops else "✗"
        print(f"  {status} {op}")
    print(f"  ✓ 算子列表加载成功")


def test_sinfon(cdo):
    """测试 4: sinfon 读取 NetCDF 文件信息"""
    separator("测试 4: sinfon — 读取文件信息")
    infiles = sorted(glob.glob(os.path.join(GPM_DATA_DIR, "GPM_Precip_*.nc")))
    if not infiles:
        print(f"  ⚠ 未找到 GPM 文件: {GPM_DATA_DIR}")
        return None

    print(f"  找到 {len(infiles)} 个 GPM 文件")
    test_file = infiles[0]
    print(f"  测试文件: {os.path.basename(test_file)}")

    info = cdo.sinfon(input=test_file)
    if isinstance(info, list):
        info = "\n".join(info)
    print(f"  sinfon 输出:\n{info[:500]}")
    print(f"  ✓ sinfon 成功")
    return infiles


def test_griddes(cdo, infiles):
    """测试 5: griddes 获取网格描述"""
    separator("测试 5: griddes — 网格描述")
    grid = cdo.griddes(input=infiles[0])
    if isinstance(grid, list):
        grid = "\n".join(grid)
    print(f"  网格信息:\n{grid[:400]}")
    print(f"  ✓ griddes 成功")


def test_remapbil(cdo, infiles):
    """测试 6: remapbil 重映射到 1°"""
    separator("测试 6: remapbil — 重映射到 1° 网格")
    outfile = os.path.join(OUTPUT_DIR, "test_remap_1deg.nc")
    if os.path.exists(outfile):
        os.remove(outfile)

    cdo.remapbil("r360x180", input=infiles[0], output=outfile)
    size_mb = os.path.getsize(outfile) / (1024 * 1024)
    print(f"  输入: {os.path.basename(infiles[0])}")
    print(f"  输出: {outfile}")
    print(f"  大小: {size_mb:.2f} MB")

    # 验证输出网格
    grid = cdo.griddes(input=outfile)
    if isinstance(grid, list):
        grid = "\n".join(grid)
    assert "xsize     = 360" in grid, "重映射网格 X 不正确"
    assert "ysize     = 180" in grid, "重映射网格 Y 不正确"
    print(f"  ✓ 重映射到 r360x180 成功 (360x180)")

    os.remove(outfile)
    return True


def test_mergetime(cdo, infiles):
    """测试 7: mergetime 合并多个时间步"""
    separator("测试 7: mergetime — 合并 3 个月")
    files_3months = infiles[:3]
    outfile = os.path.join(OUTPUT_DIR, "test_merge_3months.nc")
    if os.path.exists(outfile):
        os.remove(outfile)

    print(f"  输入文件:")
    for f in files_3months:
        print(f"    {os.path.basename(f)}")

    cdo.mergetime(input=files_3months, output=outfile)
    size_mb = os.path.getsize(outfile) / (1024 * 1024)
    print(f"  输出: {outfile}")
    print(f"  大小: {size_mb:.2f} MB")
    print(f"  ✓ mergetime 成功")

    os.remove(outfile)
    return True


def test_sellonlatbox(cdo, infiles):
    """测试 8: sellonlatbox 裁剪中国区域"""
    separator("测试 8: sellonlatbox — 裁剪中国区域")
    outfile = os.path.join(OUTPUT_DIR, "test_china_region.nc")
    if os.path.exists(outfile):
        os.remove(outfile)

    # 中国大致范围: 73-135E, 18-54N
    cdo.sellonlatbox("73,135,18,54", input=infiles[0], output=outfile)
    size_mb = os.path.getsize(outfile) / (1024 * 1024)
    print(f"  区域: 73-135°E, 18-54°N (中国)")
    print(f"  输出: {outfile}")
    print(f"  大小: {size_mb:.2f} MB")
    print(f"  ✓ sellonlatbox 裁剪成功")

    os.remove(outfile)
    return True


def test_pipeline(cdo, infiles):
    """测试 9: 管道组合 - mergetime + remapbil"""
    separator("测试 9: 管道 — mergetime + remapbil")
    files_3months = infiles[:3]
    outfile = os.path.join(OUTPUT_DIR, "test_pipeline_merge_remap.nc")
    if os.path.exists(outfile):
        os.remove(outfile)

    # CDO 支持 -op1 -op2 input 的管道链
    cdo.remapbil(
        "r360x180",
        input="-mergetime " + " ".join(files_3months),
        output=outfile
    )
    size_mb = os.path.getsize(outfile) / (1024 * 1024)
    print(f"  管道: remapbil(r360x180, mergetime(3 files))")
    print(f"  输出: {outfile}")
    print(f"  大小: {size_mb:.2f} MB")
    print(f"  ✓ 管道操作成功")

    os.remove(outfile)
    return True


def test_tempfile_output(cdo, infiles):
    """测试 10: 使用临时文件作为输出"""
    separator("测试 10: 临时文件输出")
    with tempfile.NamedTemporaryFile(suffix=".nc", delete=False) as tmp:
        tmpfile = tmp.name
    try:
        cdo.timmean(input=infiles[0], output=tmpfile)
        size_kb = os.path.getsize(tmpfile) / 1024
        print(f"  timmean 输出: {tmpfile}")
        print(f"  临时文件大小: {size_kb:.1f} KB")
        print(f"  ✓ 临时文件输出成功")
    finally:
        if os.path.exists(tmpfile):
            os.remove(tmpfile)


def main():
    print("=" * 60)
    print("  skyborn-cdo 综合测试")
    print(f"  Python {sys.version}")
    print("=" * 60)

    passed = 0
    failed = 0
    total = 0

    try:
        # 基础测试（无需数据文件）
        cdo = test_import()
        total += 1
        passed += 1

        test_version(cdo)
        total += 1
        passed += 1

        test_operators(cdo)
        total += 1
        passed += 1

        # 数据测试（需要 GPM 文件）
        infiles = test_sinfon(cdo)
        total += 1
        passed += 1

        if infiles:
            test_griddes(cdo, infiles)
            total += 1
            passed += 1

            test_remapbil(cdo, infiles)
            total += 1
            passed += 1

            test_mergetime(cdo, infiles)
            total += 1
            passed += 1

            test_sellonlatbox(cdo, infiles)
            total += 1
            passed += 1

            test_pipeline(cdo, infiles)
            total += 1
            passed += 1

            test_tempfile_output(cdo, infiles)
            total += 1
            passed += 1

    except Exception as e:
        failed += 1
        total += 1
        print(f"\n  ✗ 测试失败: {e}")
        import traceback
        traceback.print_exc()

    separator("测试结果")
    print(f"  通过: {passed}/{total}")
    if failed:
        print(f"  失败: {failed}/{total}")
    else:
        print(f"  全部通过! ✓")

    return 0 if failed == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
