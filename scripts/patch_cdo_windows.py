#!/usr/bin/env python3
"""
Smart Windows Compatibility Patcher for CDO Source Code

使用代码特征匹配（而非固定行号）自动修改 CDO 源码，适配 Windows 编译环境。
CDO 版本更新后，此脚本大概率仍能自动工作，降低维护成本。

用法:
    python patch_cdo_windows.py apply [--cdo-src PATH]    # 应用修改
    python patch_cdo_windows.py restore [--cdo-src PATH]  # 恢复原始代码
    python patch_cdo_windows.py verify [--cdo-src PATH]   # 验证（不修改）
"""

import re
import sys
import argparse
from pathlib import Path
from typing import Tuple


class WindowsPatcher:
    """CDO Windows 兼容性补丁管理器"""
    
    def __init__(self, cdo_src_dir: Path):
        self.cdo_src = Path(cdo_src_dir).resolve()
        self.backup_dir = self.cdo_src / ".patch_backup"
        
    def patch_file(self, rel_path: str, patches: list, dry_run: bool = False) -> Tuple[bool, int]:
        """对单个文件应用补丁。返回 (success, count)"""
        file_path = self.cdo_src / rel_path
        
        if not file_path.exists():
            print(f"[X] {rel_path}: 文件不存在")
            return False, 0
        
        print(f"[*] {rel_path}")
        
        try:
            content = file_path.read_text(encoding='utf-8', errors='ignore')
        except Exception as e:
            print(f"   [X] 读取失败: {e}")
            return False, 0
        
        original = content
        applied_count = 0
        
        for desc, pattern, replacement in patches:
            if isinstance(pattern, str):
                # 简单字符串替换
                if pattern in content:
                    content = content.replace(pattern, replacement, 1)
                    print(f"   [+] {desc}")
                    applied_count += 1
                else:
                    print(f"   [ ] 未找到: {desc}")
            else:
                # 正则表达式替换
                new_content, count = pattern.subn(replacement, content, count=1)
                if count > 0:
                    content = new_content
                    print(f"   [+] {desc}")
                    applied_count += 1
                else:
                    print(f"   [ ] 未找到: {desc}")
        
        # 写入修改
        if content != original and not dry_run:
            # 备份原文件
            backup_path = self.backup_dir / rel_path.replace('/', '_').replace('\\', '_')
            backup_path.parent.mkdir(parents=True, exist_ok=True)
            backup_path.write_text(original, encoding='utf-8', newline='\n')
            
            # 写入修改后的内容
            file_path.write_text(content, encoding='utf-8', newline='\n')
        
        return applied_count > 0, applied_count
    
    def apply_all(self, dry_run: bool = False) -> bool:
        """应用所有补丁"""
        print(f"CDO 源码目录: {self.cdo_src}\n")
        
        if not dry_run:
            self.backup_dir.mkdir(exist_ok=True)
        
        total_files = 0
        total_patches = 0
        
        # =================================================================
        # 补丁定义：基于代码特征匹配，不依赖固定行号
        # =================================================================
        
        patches = [
            # --- src/cdo.cc ---
            ("src/cdo.cc", [
                ("添加 Windows 头文件 (io.h, windows.h)",
                 re.compile(r'^(#include\s+<unistd\.h>.*?)$', re.MULTILINE),
                 r'#ifdef _WIN32\n#include <io.h>\n#include <windows.h>\n#else\n\1\n#endif'),
                
                ("cdo_init_is_tty: Windows 实现",
                 re.compile(
                     r'(static\s+void\s+cdo_init_is_tty\s*\(\s*\)\s*\{)'
                     r'([^}]+)'
                     r'(\})',
                     re.DOTALL
                 ),
                 lambda m: (
                     m.group(1) + '\n#ifdef _WIN32\n'
                     '  cdo::stdinIsTerminal = _isatty(_fileno(stdin));\n'
                     '  cdo::stdoutIsTerminal = _isatty(_fileno(stdout));\n'
                     '  cdo::stderrIsTerminal = _isatty(_fileno(stderr));\n'
                     '#else' + m.group(2) + '#endif\n' + m.group(3)
                 )),
                
                ("在 clear_processes 前添加 fflush",
                 re.compile(r'(\s+)(g_processManager\.clear_processes\s*\(\s*\)\s*;)'),
                 r'\1fflush(stdout);\n\1\2'),
            ]),
            
            # --- src/cdo_getopt.cc ---
            ("src/cdo_getopt.cc", [
                ("屏蔽 sys/ioctl.h (Windows 不可用)",
                 re.compile(r'^(#include\s+<sys/ioctl\.h>)$', re.MULTILINE),
                 r'#ifndef _WIN32\n\1\n#endif'),
                
                ("Screening unistd.h", re.compile(r'^(#include\s+<unistd\.h>)$', re.MULTILINE),
                 r'#ifndef _WIN32\n\1\n#endif'),
            ]),
            
            # --- 其他文件：统一屏蔽 unistd.h ---
            *[(f, [("屏蔽 unistd.h",
                   re.compile(r'^(#include\s+<unistd\.h>)$', re.MULTILINE),
                   r'#ifndef _WIN32\n\1\n#endif')])
              for f in [
                  "src/cdo_zaxis.cc",
                  "src/dcw_reader.cc",
                  "src/expr_lex.cc",
                  "src/griddes.cc",
                  "src/merge_axis.cc",
                  "src/operators/CMOR.cc",
              ]],
        ]
        
        # 应用所有补丁
        for rel_path, file_patches in patches:
            success, count = self.patch_file(rel_path, file_patches, dry_run)
            if success:
                total_files += 1
                total_patches += count
            print()
        
        # 总结
        mode = "验证模式 - 未修改文件" if dry_run else f"{total_files} 个文件已修改"
        print(f"{'='*70}")
        print(f"完成: {mode}, 共应用 {total_patches} 处补丁")
        print(f"{'='*70}")
        
        return total_files > 0
    
    def restore_all(self) -> bool:
        """从备份恢复所有文件"""
        if not self.backup_dir.exists():
            print("未找到备份 - 无需恢复")
            return True
        
        print(f"从备份恢复: {self.backup_dir}\n")
        
        restored = 0
        for backup_file in self.backup_dir.iterdir():
            if backup_file.is_file():
                # 转换备份文件名回原路径
                rel_path = backup_file.name.replace('_', '/', 2)  # src_cdo.cc -> src/cdo.cc
                original_file = self.cdo_src / rel_path
                
                try:
                    content = backup_file.read_text(encoding='utf-8')
                    original_file.write_text(content, encoding='utf-8', newline='\n')
                    print(f"[+] 已恢复: {rel_path}")
                    restored += 1
                except Exception as e:
                    print(f"[X] 恢复失败 {rel_path}: {e}")
        
        # 清理备份目录
        try:
            for f in self.backup_dir.iterdir():
                f.unlink()
            self.backup_dir.rmdir()
            print(f"\n[OK] {restored} 个文件已恢复，备份目录已删除")
        except Exception as e:
            print(f"\n[!] {restored} 个文件已恢复，但无法删除备份目录: {e}")
        
        return True


def main():
    parser = argparse.ArgumentParser(description="CDO Windows 智能补丁工具")
    parser.add_argument("action", choices=["apply", "restore", "verify"], 
                        help="操作: apply=应用补丁, restore=恢复原始, verify=验证")
    parser.add_argument("--cdo-src", type=Path,
                        help="CDO 源码目录 (默认: ../vendor/cdo)")
    
    args = parser.parse_args()
    
    # 确定源码目录
    if args.cdo_src:
        cdo_src = args.cdo_src
    else:
        script_dir = Path(__file__).parent
        cdo_src = script_dir.parent / "vendor" / "cdo"
    
    if not cdo_src.exists():
        print(f"[X] CDO 源码目录不存在: {cdo_src}")
        print("   请使用 --cdo-src 指定或确保 vendor/cdo 存在")
        return 1
    
    patcher = WindowsPatcher(cdo_src)
    
    if args.action == "apply":
        success = patcher.apply_all(dry_run=False)
    elif args.action == "restore":
        success = patcher.restore_all()
    elif args.action == "verify":
        success = patcher.apply_all(dry_run=True)
    
    return 0 if success else 1


if __name__ == "__main__":
    sys.exit(main())
