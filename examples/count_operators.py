"""
Query CDO Operator Count and Classification

Demonstrates how to use skyborn-cdo to query the available 938 CDO operators
and group them by category.
"""
from collections import defaultdict
from skyborn_cdo import Cdo

cdo = Cdo()

# Get all operators
operators = cdo.operators()
print(f"Total CDO Operators: {len(operators)}")

# Group and count by first letter
by_prefix = defaultdict(list)
for op in sorted(operators):
    by_prefix[op[0]].append(op)

print("\nDistribution by First Letter:")
for letter in sorted(by_prefix.keys()):
    print(f"  {letter}: {len(by_prefix[letter])} operators")

# List operators in commonly used categories
print("\n" + "=" * 70)
print("Common Operator Categories")
print("=" * 70)

# Time operations
time_ops = [op for op in operators if any(
    x in op for x in ['time', 'mon', 'year', 'day', 'hour', 'seas'])]
print(f"\nTime Operations ({len(time_ops)} operators):")
print(f"  Examples: {', '.join(sorted(time_ops)[:12])}")

# Spatial selection
sel_ops = [op for op in operators if op.startswith('sel')]
print(f"\nSpatial Selection ({len(sel_ops)} operators):")
print(f"  {', '.join(sorted(sel_ops)[:15])}")

# Statistics
stat_ops = [op for op in operators if any(
    x in op for x in ['mean', 'std', 'var', 'min', 'max', 'sum'])]
print(f"\nStatistics Operators ({len(stat_ops)} operators):")
print(f"  Examples: {', '.join(sorted(stat_ops)[:15])}")

# Remapping
remap_ops = [op for op in operators if 'remap' in op]
print(f"\nRemapping ({len(remap_ops)} operators):")
print(f"  {', '.join(sorted(remap_ops))}")

# Arithmetic
arith_ops = [op for op in operators if any(
    x in op for x in ['add', 'sub', 'mul', 'div', 'sqrt', 'exp', 'log', 'abs'])]
print(f"\nArithmetic Operations ({len(arith_ops)} operators):")
print(f"  Examples: {', '.join(sorted(arith_ops)[:20])}")

# Format conversion
format_ops = [op for op in operators if any(
    x in op for x in ['copy', 'import', 'output', 'split'])]
print(f"\nFormat/Conversion ({len(format_ops)} operators):")
print(f"  Examples: {', '.join(sorted(format_ops)[:15])}")

# Info class
info_ops = [op for op in operators if any(
    x in op for x in ['show', 'info', 'print', 'ncode', 'ntime', 'nlevel'])]
print(f"\nInfo/Query ({len(info_ops)} operators):")
print(f"  Examples: {', '.join(sorted(info_ops)[:20])}")

# Vertical layer operations
vert_ops = [op for op in operators if any(
    x in op for x in ['ml2', 'ml2pl', 'intlevel', 'pressure'])]
print(f"\nVertical Layer Operations ({len(vert_ops)} operators):")
print(f"  {', '.join(sorted(vert_ops)[:15])}")

# Grid operations
grid_ops = [op for op in operators if any(
    x in op for x in ['grid', 'setgrid', 'gp2', 'sp2'])]
print(f"\nGrid Operations ({len(grid_ops)} operators):")
print(f"  Examples: {', '.join(sorted(grid_ops)[:15])}")

print("\n" + "=" * 70)
print("View Help for Specific Operators")
print("=" * 70)
print("\nUsage:")
print("  cdo.help('operator_name')  # View detailed help for any operator")
print("  cdo.has_operator('xxx')     # Check if operator exists")
print("\nExample:")
print("\n>>> cdo.help('sellonlatbox')")
print(cdo.help('sellonlatbox')[:400] + "...")
