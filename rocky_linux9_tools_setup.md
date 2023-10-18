## Installing all the tools required in Rocky Linux 9
### Install Rocky Linux 9
Download and install [Rocky 9](https://rockylinux.org/download/) minimal in a server. [VMware](https://www.vmware.com/) virtual machine will also be feesible.

### Dependancies
```bash
yum groupinstall "Development Tools"
yum install -y openssl-devel libuuid-devel libseccomp-devel wget squashfs-tools go go-devel python-pip bzip2 bzip2-devel ncurses-devel curl libcurl libcurl-devel cpan
cpan install Sys::Hostname
pip install ninja
```

### Singularity
```bash
wget https://github.com/sylabs/singularity/releases/download/v4.0.1/singularity-ce-4.0.1-1.el9.x86_64.rpm
yum localinstall singularity-ce-4.0.1-1.el9.x86_64.rpm
```

### htslib
```bash
wget https://github.com/samtools/htslib/releases/download/1.18/htslib-1.18.tar.bz2
tar -xvf htslib-1.18.tar.bz2 
cd htslib-1.18/
./configure
make
make install
export LD_LIBRARY_PATH=/usr/local/lib
```

### SAMtools
```bash
wget https://github.com/samtools/samtools/releases/download/1.18/samtools-1.18.tar.bz2
tar -xvf samtools-1.18.tar.bz2 
cd samtools-1.18/
./configure
make
make install
```

### BCFtools
```bash
wget https://github.com/samtools/bcftools/releases/download/1.18/bcftools-1.18.tar.bz2
tar -xvf bcftools-1.18.tar.bz2 
cd bcftools-1.18/
./configure
make
make install
```

### Bowtie-2
```bash
wget https://github.com/BenLangmead/bowtie2/archive/refs/tags/v2.5.2.tar.gz
tar -xvf v2.5.2.tar.gz 
cd bowtie2-2.5.2
make
make install
```

### VG
```bash
wget https://github.com/vgteam/vg/releases/download/v1.51.0/vg
chmod 755 vg
cp vg /usr/local/bin/
```

### Freebayes
```bash
git clone --recursive https://github.com/freebayes/freebayes.git
cd freebayes/
meson build/ --buildtype debug
cd build/
ninja
cp freebayes /usr/local/bin/
```

### NGSNGS
```bash
git clone https://github.com/RAHenriksen/NGSNGS.git
cd NGSNGS
make
./ngsngs 
cp ngsngs /usr/local/bin/ngsngs
```

### BWA
```bash
wget https://github.com/lh3/bwa/releases/download/v0.7.17/bwa-0.7.17.tar.bz2
tar -xvf bwa-0.7.17.tar.bz2 
cd bwa-0.7.17/
vim rle.h # Then comment the line 33 -> // const uint8_t rle_auxtab[8];
make
cp bwa /usr/local/bin/
```

### Singularity Images
```bash
singularity pull docker://broadinstitute/gatk:4.4.0.0
singularity pull docker://google/deepvariant:1.5.0
```
