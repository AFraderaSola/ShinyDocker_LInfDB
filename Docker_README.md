##### Dockerize a Shiny App ##########

I mainly followed the instructions on the following blog:

https://www.r-bloggers.com/2021/05/dockerizing-shiny-applications/

If you want to learn more, just read there.

Below I will just detail the steps I used to dockerize the Leishmania infectome App.

### Summary ###

1 - Adapt the Dockerfile document
2 - Run the Dockerfile document to create a docker image
3 - Push the image to a repository
4 - Pull the image on any other system, the app should be ready to run

### 1 - Adapt the Dockerfile document ###

I work on the following directory:

Home/Desktop/ShinyToDocker/LeishmaniaInfectome

There, I have the following files and directorys:

1- "app" directory, containing all the shiny app files

2- "RLibraries" directory, containing all libraries that we want to install manually, instead from the cran/bioconductor repository (e.g., rMQanalysis library)

2- "Dockerfile" document, containing the instructions to dockerize the app directory

Here are some of comments on how to adapt the docker file document to each app:

##

# With this command, you detail which r-base version you want #

FROM rocker/r-base:latest

##

LABEL maintainer="AFraderaSola <A.FraderaSola@imb-mainz.de>"

##

# With this command, you install any ubuntu libraries needed to run R and their packages.#

RUN apt-get update && apt-get install -y --no-install-recommends \
    sudo \
    default-jdk \
    r-cran-rjava \
    libcurl4-gnutls-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libxml2-dev \
    libcairo2-dev \
    libxt-dev \
    libssl-dev \
    cmake \
    libssh2-1-dev \
    libpng-dev \ 
    libtiff5-dev \ 
    libjpeg-dev \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*
    
##  

# With this command, you install the R libraries you need for the shiny app #
# You need to specify the working directory with WORKDIR #
  
RUN R -e "install.packages(c('shiny', 'ggplot2', 'tidyverse', 'ggpubr', 'plotly', 'ggrepel', 'shinyBS', 'DT', 'mailtoR', 'network', 'scales', 'sna', 'intergraph', 'ggparty', 'ggraph', 'igraph', 'visNetwork', 'shinythemes','nloptr', 'kohonen', 'effectsize', 'gplots', 'Peptides', 'maptools', 'seqinr', 'xlsx', 'logspline', 'splitstackshape'), dependencies=TRUE, repos='http://cran.rstudio.com/')"
WORKDIR /home/Desktop/ShinyToDocker/LeishmaniaInfectome/
COPY ./RLibraries/cfpscripts_0.4.11.9008.tar.gz .
COPY ./RLibraries/rMQanalysis_0.3.4.9012.tar.gz .
RUN R -e "install.packages(c('rMQanalysis_0.3.4.9012.tar.gz', 'cfpscripts_0.4.11.9008.tar.gz'),repos = NULL, type = 'source')"

##

RUN echo "local(options(shiny.port = 3838, shiny.host = '0.0.0.0'))" > /usr/lib/R/etc/Rprofile.site
RUN addgroup --system app \
    && adduser --system --ingroup app app
    
##  

# With this command, you change the working directory to the app directory #
    
WORKDIR /home/Desktop/ShinyToDocker/LeishmaniaInfectome/app
COPY app .
RUN chown app:app -R /home/Desktop/ShinyToDocker/LeishmaniaInfectome/app
USER app

##  

# I have not tried this, but you can probably change here which port to send the app to #

EXPOSE 3838

##

CMD ["R", "-e", "shiny::runApp('/home/Desktop/ShinyToDocker/LeishmaniaInfectome/app')"]/app')"]


### 2 - Run the Dockerfile document ###

# You can build the docker image with the following command: #

sudo docker build -t ghcr.io/afraderasola/shinydocker/leishmania_infectome .

# ghcr.io/afraderasola/ takes to my github page; there are other options where to allocate the image, like gitlab. You then should check that the image build work but running it with the following command: #

sudo docker run -p 3838:3838 ghcr.io/afraderasola/shinydocker/leishmania_infectome

# If the image works, you should get a "listening" message and you should be able to open the app on a browser #

### 3 - Push the image to a repository ###

# If the image works as expected, you can push it the repository with the following command: #

sudo docker push ghcr.io/afraderasola/shinydocker/leishmania_infectome:latest

### 4 - Pull the image to another system ###

# Once the image is on the repository, you can pull it to another system (e.g., the butterVM server) and get it running there:

sudo docker pull ghcr.io/afraderasola/shinydocker/leishmania_infectome:latest

# That's it! You should now have your new docker shiny app running! #
