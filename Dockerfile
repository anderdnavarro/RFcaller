FROM debian:buster-slim
## Install linux packages and remove cache
RUN ["/bin/bash", "-c", "apt-get update && apt-get install -y apt-utils && apt-get install -y bc bzip2 curl dialog gawk gcc git jq less libbz2-dev libcairo2-dev libcurl4-openssl-dev libgif-dev libjpeg-dev liblzma-dev libncurses5-dev libncursesw5-dev libssl-dev make nano procps python3 python3-pip unzip wget zip zlib1g-dev && apt-get autoremove -y && apt-get clean && rm -rf /var/lib/apt/lists/*"]
## Install python modules
RUN ["/bin/bash", "-c", "pip3 install --no-cache-dir argparse==1.4.0 asn1crypto==0.24.0 certifi==2020.6.20 chardet==3.0.4 click==7.1.2 cloudpickle==1.5.0 cryptography==2.6.1 cycler==0.10.0 Cython==0.29.17 dask-glm==0.2.0 dask-ml==1.3.0 dask==2.15.0 deap==1.3.1 distributed==2.22.0 entrypoints==0.3 fsspec==0.7.3 HeapDict==1.0.1 idna==2.10 joblib==0.14.1 keyring==17.1.1 keyrings.alt==3.1.1 kiwisolver==1.2.0 llvmlite==0.33.0 locket==0.2.0 matplotlib==3.3.0 msgpack==1.0.0 multipledispatch==0.6.0 numba==0.50.1 numpy==1.18.4 packaging==20.4 pandas==1.0.3 partd==1.1.0 Pillow==7.2.0 psutil==5.7.2 pycairo==1.19.1 pycrypto==2.6.1 PyGObject==3.30.4 pyparsing==2.4.7 pytabix==0.1 python-dateutil==2.8.1 pytz==2020.1 pyxdg==0.25 PyYAML==5.3.1 regex==2020.4.4 requests==2.24.0 scikit-learn==0.22.2.post1 scikit-MDR==0.4.4 scipy==1.4.1 SecretStorage==2.3.1 six==1.12.0 sklearn==0.0 skrebate==0.6 sortedcontainers==2.2.2 stopit==1.1.2 tblib==1.7.0 toolz==0.10.0 tornado==6.0.4 TPOT==0.11.1 tqdm==4.46.0 update-checker==0.17 urllib3==1.25.10 xgboost==1.0.2 zict==2.0.0 && pip3 install --no-cache-dir pysam==0.15.4"]
## Install samtools
COPY samtools/ /home/samtools/
RUN ["/bin/bash", "-c", "cd /home/samtools; for i in htslib-1.10.2 bcftools-1.10.2 samtools-1.10; do tar xvjf ${i}.tar.bz2; cd ${i}; ./configure; make; make install; cd ../; done; cd /home/; rm -rf samtools"]
## Install bedtools
COPY bedtools.zip /usr/local/bin/
RUN ["/bin/bash", "-c", "unzip -qq bedtools.zip; rm bedtools.zip"]
## Copy necessary files
COPY databases/ /home/databases/
COPY RFcaller/ /home/RFcaller/
RUN ["/bin/bash", "-c", "mkdir /output"]
## Entrypoint
ENTRYPOINT ["RFcaller"]
## Set environment
ENV TZ=Europe/Madrid
ENV PATH=$PATH:/home/RFcaller/scripts/
WORKDIR /output
## Labels
LABEL Author="Ander Díaz-Navarro, Pablo Bousquets-Muñoz & Xose S. Puente - Universidad de Oviedo"
LABEL Contact="ander.d.navarro@gmail.com"
