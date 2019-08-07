FROM ubuntu

RUN apt-get update && apt-get install -y \
    apt-utils make 
    
RUN apt-get install -y  --fix-missing \
    gfortran xorg-dev 

RUN apt-get install -y  --fix-missing \
    libreadline-dev 

RUN apt-get install -y --fix-missing \
    libxft2
    
RUN apt-get install -y --fix-missing \
    groff 
    
RUN apt-get install -y --fix-missing \
    x11-apps

RUN apt-get install -y --fix-missing \
    imagemagick

RUN apt-get update && apt-get install -y  pdf2svg

RUN apt-get update && apt-get install -y  texlive 

RUN apt-get update && apt-get install -y  texlive-extra-utils

RUN apt-get update && apt-get install -y  pandoc

RUN apt-get update && apt-get install -y  python-bs4

RUN rm -rf /var/lib/apt/lists/*

WORKDIR /opt/SAD/SRC

COPY . .

RUN make mostlyclean
    
RUN make install clean

WORKDIR ../

RUN rm -r ./SRC

CMD ["/opt/SAD/bin/gs","/opt/SAD/examples/design_example.sad"]
