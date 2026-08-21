import sys
import importlib
import subprocess
import shutil
import os

def get_tool_version(command, version_flag='--version'):
    try:
        executable = command.split()[0]
        if executable != 'java' and shutil.which(executable) is None and not os.path.exists(executable):
             return "not_found"
        
        try:
            output = subprocess.check_output(f"{command} {version_flag}", shell=True, stderr=subprocess.STDOUT).decode().strip()
        except subprocess.CalledProcessError as e:
            output = e.output.decode().strip()
            
        lines = output.split('\n')
        
        if 'ivar' in command:
            for line in lines:
                if 'iVar version' in line:
                    return line.split('iVar version')[1].strip()

        for line in lines:
            if "Looking to launch" in line:
                continue
            if any(c.isdigit() for c in line):
                return line.strip().replace('"', '')
                
        return lines[0].strip()
    except Exception as e:
        return f"error"

def main():
    packages = [
        "pandas", "numpy", "pysam", "Bio"
    ]

    print("python_tools:")
    print(f"  python: {sys.version.split()[0]}")

    for p in packages:
        try:
            mod = importlib.import_module(p)
            version = getattr(mod, "__version__", "unknown")
            print(f"  {p}: {version}")
        except ImportError:
            print(f"  {p}: not_installed")

    print("system_tools:")

    tools = {
        "bcftools": ("bcftools", "--version"),
        "samtools": ("samtools", "--version"),
        "bowtie2": ("bowtie2", "--version"),
        "minimap2": ("minimap2", "--version"),
        "fastp": ("fastp", "--version"),
        "fastqc": ("fastqc", "--version"),
        "multiqc": ("multiqc", "--version"),
        "chopper": ("chopper", "--version"),
        "medaka": ("medaka", "--version"),
        "ivar": ("ivar", "version"),
        "bedtools": ("bedtools", "--version"),
        "mosdepth": ("mosdepth", "--version"),
        "snpeff": ("snpEff", "-version"),
        "snpsift": ("snpSift", ""),
        "nextclade": ("nextclade", "--version"),
        "hostile": ("hostile", "--version"),
        "megahit": ("megahit", "--version"),
        "cutadapt": ("cutadapt", "--version"),
        "quast": ("quast.py", "--version"),
        "blastn": ("blastn", "-version"),
        "picard": ("PicardCommandLine", "MarkDuplicates --version"),
        "java": ("java", "-version"),
    }

    for name, (cmd, flag) in tools.items():
        ver = get_tool_version(cmd, flag)
        print(f"  {name}: \"{ver}\"")

if __name__ == "__main__":
    main()