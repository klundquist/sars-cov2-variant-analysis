# Use an official Ubuntu as a parent image
FROM ubuntu:20.04

# Set environment variables
ENV DEBIAN_FRONTEND=noninteractive

# Install essential packages
RUN apt-get update && apt-get install -y \
    curl \
    vim \
    unzip \
    python3 \
    python3-pip \
    wget \
    tar \
    cmake \
    g++ \
    && apt-get clean

# Install Miniconda for ARM64 architecture if available or x86_64 otherwise
RUN wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O /tmp/miniconda.sh && \    
    bash /tmp/miniconda.sh -b -p /opt/conda && \
    rm /tmp/miniconda.sh && \
    /opt/conda/bin/conda clean -a -y && \
    ln -s /opt/conda/bin/conda /usr/local/bin/conda

# Download and extract the SRA Toolkit
RUN wget -O /tmp/sratoolkit.tar.gz https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/current/sratoolkit.current-ubuntu64.tar.gz && \
    mkdir -p /opt/sratoolkit && \
    tar -vxzf /tmp/sratoolkit.tar.gz -C /opt/sratoolkit --strip-components=1 

# Add the SRA Toolkit binary directory to the PATH
ENV PATH="/opt/sratoolkit/bin:$PATH"

# Configure conda and install bioinformatics tools
RUN conda config --add channels defaults && \
    conda config --add channels bioconda && \
    conda config --add channels conda-forge && \
    conda install -y \
    bwa \
    samtools \
    vcftools \
    fastqc

# Install Java
RUN apt-get update && apt-get install -y \
    openjdk-17-jdk

# Install GATK
RUN wget https://github.com/broadinstitute/gatk/releases/download/4.6.0.0/gatk-4.6.0.0.zip && \
    unzip gatk-4.6.0.0.zip && \
    mv gatk-4.6.0.0 /opt/gatk && \
    ln -s /opt/gatk/gatk /usr/bin/gatk

# Install SnpEff
RUN wget https://snpeff.blob.core.windows.net/versions/snpEff_latest_core.zip && \
    unzip snpEff_latest_core.zip -d /opt/ && \
    echo "alias snpEff='java -jar /opt/snpEff/snpEff.jar'" >> ~/.bash_aliases

# Install Trimmomatic
RUN wget http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.39.zip && \
    unzip Trimmomatic-0.39.zip -d /opt/ && \
    echo "alias trimmomatic='java -jar /opt/Trimmomatic-0.39/trimmomatic-0.39.jar'" >> ~/.bash_aliases

# Install FreeBayes manually
RUN wget https://github.com/freebayes/freebayes/releases/download/v1.3.6/freebayes-1.3.6-linux-amd64-static.gz && \
    gunzip freebayes-1.3.6-linux-amd64-static.gz && \
    chmod +x freebayes-1.3.6-linux-amd64-static && \
    mv freebayes-1.3.6-linux-amd64-static /usr/local/bin/freebayes

# Install MultiQC
RUN pip3 install multiqc

# Make sure that .bash_aliases is sourced by .bashrc
RUN echo "source ~/.bash_aliases" >> ~/.bashrc

# Set the working directory
WORKDIR /data

# Set PATH
ENV PATH="/opt/gatk:/opt/snpEff:/opt/conda/bin:${PATH}"
# Set JAVA_HOME and update PATH
ENV JAVA_HOME=/usr/lib/jvm/java-17-openjdk-amd64
ENV PATH="$JAVA_HOME/bin:$PATH"

# Define default command
CMD ["/bin/bash"]
