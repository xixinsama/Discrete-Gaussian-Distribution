from elftools.elf.elffile import ELFFile

def get_elf_section_sizes(elf_file_path):
    with open(elf_file_path, 'rb') as f:
        elf = ELFFile(f)
        sizes = {
            '.text': 0,
            '.data': 0,
            '.bss': 0,
            '.rodata': 0
        }
        for section in elf.iter_sections():
            if section.name in sizes:
                sizes[section.name] = section['sh_size']
        return sizes

# 替换 'your_program.elf' 为你的可执行文件路径
elf_sizes = get_elf_section_sizes('/home/luoluo/crypto/saiti2/FinalA/example')
for section, size in elf_sizes.items():
    print(f'{section} size: {size} bytes')
