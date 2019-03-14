
##############################################################
# Dockerfile Version:   1.8
# Software:             STARK-BASE
# Software Version:     0.9.18b
# Software Website:     none
# Licence:              GNU Affero General Public License (AGPL)
# Description:          STARK
# Usage:                docker run -ti [-v [DATA FOLDER]:/data -v [DATABASE_FOLDER]:/databases] stark-base:version
##############################################################

##########
# README #
##########

# Config parameters
#    identify yum packages for installation
#    identify yum packages to remove
#
# Dependecies installation
#    identify tools dependences
#    config each tool
#    write installation procedure for each tools
#
# Tool
#    configure tool
#    write isntallation procedure for the tool
#    add link to current and root tool folder
#
# Workdir / Entrypoint / Cmd
#    configure workdir, endpoint and command
#    /!\ no variables in endpoint



########
# FROM #
########

FROM centos:7
MAINTAINER Antony Le Bechec <antony.lebechec@gmail.com>
LABEL Software="STARK-BASE" \
	Version="1.3" \
	Website="none" \
	Description="STARK-BASE" \
	License="GNU Affero General Public License (AGPL)" \
	Usage="docker run -ti [-v [DATA FOLDER]:/STARK/data -v [DATABASE_FOLDER]:/STARK/databases] stark:version"



##############
# PARAMETERS #
##############

ENV STARK_FOLDER=/STARK
ENV TOOLS=$STARK_FOLDER/tools
ENV DATA=$STARK_FOLDER/data
ENV TOOL=$STARK_FOLDER/tool
ENV DATABASES=$STARK_FOLDER/databases
ENV YUM_INSTALL="autoconf automake bc bzip2 bzip2-devel curl gcc gcc-c++ git java java-1.8.0 lzma lzma-devel make ncurses-devel perl perl-Data-Dumper perl-Digest-MD5 perl-Switch perl-devel perl-Tk tbb-devel unzip wget which xz xz-devel zlib zlib-devel zlib2 zlib2-devel ghostscript enscript"
ENV YUM_REMOVE="autoconf automake bzip2-devel lzma-devel ncurses-devel perl-devel tbb-devel xz-devel zlib-devel zlib2-devel"

#epel-release R


###############
# YUM INSTALL #
###############

RUN yum update -y ; \
	yum install -y epel-release ; \
	yum install -y $YUM_INSTALL ;



################
# DEPENDENCIES #
################


###########
# ANNOVAR #
###########

ENV TOOL_NAME=annovar
ENV TOOL_VERSION=2018Apr16
ENV TARBALL_LOCATION=http://www.openbioinformatics.org/annovar/download/0wgxR2rIVP/
ENV TARBALL=annovar.latest.tar.gz
ENV TARBALL_FOLDER=$TOOL_NAME
ENV TOOL_DATABASE_FOLDER=$TOOLS/$TOOL_NAME/$TOOL_VERSION/bin/databases/
ENV TOOL_DATABASE_FOLDER_LINK=$DATABASES/annovar_sources
ENV DEST=$TOOLS/$TOOL_NAME/$TOOL_VERSION
ENV PATH=$TOOLS/$TOOL_NAME/$TOOL_VERSION/bin:$PATH

RUN wget $TARBALL_LOCATION/$TARBALL ; \
    tar xf $TARBALL ; \
    rm -rf $TARBALL ; \
    cd $TARBALL_FOLDER ; \
    mkdir -p $TOOLS/$TOOL_NAME/$TOOL_VERSION/bin ; \
    cp *.pl $TOOLS/$TOOL_NAME/$TOOL_VERSION/bin/ -R ; \
    cd ../ ; \
    rm -rf $TARBALL_FOLDER ; \
    ln -s $TOOLS/$TOOL_NAME/$TOOL_VERSION $TOOLS/$TOOL_NAME/current ; \
	mkdir -p $TOOL_DATABASE_FOLDER_LINK ; \
	mkdir -p $TOOL_DATABASE_FOLDER ; \
	ln -s $TOOL_DATABASE_FOLDER_LINK $TOOL_DATABASE_FOLDER ;


##########
# HTSLIB #
##########

ENV TOOL_NAME=htslib
ENV TOOL_VERSION=1.8
ENV TARBALL_LOCATION=https://github.com/samtools/$TOOL_NAME/releases/download/$TOOL_VERSION/
ENV TARBALL=$TOOL_NAME-$TOOL_VERSION.tar.bz2
ENV DEST=$TOOLS/$TOOL_NAME/$TOOL_VERSION
ENV PATH=$TOOLS/$TOOL_NAME/$TOOL_VERSION/bin:$PATH

# INSTALL
RUN wget $TARBALL_LOCATION/$TARBALL ; \
    tar xf $TARBALL ; \
    rm -rf $TARBALL ; \
    cd $TOOL_NAME-$TOOL_VERSION ; \
    make prefix=$TOOLS/$TOOL_NAME/$TOOL_VERSION install ; \
	ln -s $TOOLS/$TOOL_NAME/$TOOL_VERSION $TOOLS/$TOOL_NAME/current ; \
    cd ../ ; \
    rm -rf $TOOL_NAME-$TOOL_VERSION ;



############
# BCFTOOLS #
############

ENV TOOL_NAME=bcftools
ENV TOOL_VERSION=1.8
ENV TARBALL_LOCATION=https://github.com/samtools/$TOOL_NAME/releases/download/$TOOL_VERSION/
ENV TARBALL=$TOOL_NAME-$TOOL_VERSION.tar.bz2
ENV DEST=$TOOLS/$TOOL_NAME/$TOOL_VERSION
ENV PATH=$TOOLS/$TOOL_NAME/$TOOL_VERSION/bin:$PATH

# INSTALL
RUN wget $TARBALL_LOCATION/$TARBALL ; \
    tar xf $TARBALL ; \
    rm -rf $TARBALL ; \
    cd $TOOL_NAME-$TOOL_VERSION ; \
    make prefix=$TOOLS/$TOOL_NAME/$TOOL_VERSION install ; \
	ln -s $TOOLS/$TOOL_NAME/$TOOL_VERSION $TOOLS/$TOOL_NAME/current ; \
    cd ../ ; \
    rm -rf $TOOL_NAME-$TOOL_VERSION ;



#############
# BCL2FASTQ #
#############

ENV TOOL_NAME=bcl2fastq
ENV TOOL_VERSION=2.20.0
ENV TOOL_VERSION_FOR_FILE=2-20-0
ENV ZIPBALL_LOCATION=https://support.illumina.com/content/dam/illumina-support/documents/downloads/software/bcl2fastq/
ENV ZIPBALL=$TOOL_NAME"2-v"$TOOL_VERSION_FOR_FILE"-linux-x86-64.zip"
ENV RPM=$TOOL_NAME"2-v"$TOOL_VERSION"*.rpm"
ENV DEST=$TOOLS/$TOOL_NAME/$TOOL_VERSION
ENV PATH=$TOOLS/$TOOL_NAME/$TOOL_VERSION/bin:$PATH
ENV TOOLS_DATA_RUN=/data/run
ENV TOOLS_DATA_OUTPUT=/data/output

RUN wget $ZIPBALL_LOCATION/$ZIPBALL ; \
    unzip  $ZIPBALL ; \
    yum -y --nogpgcheck localinstall $RPM ; \
    rm -f $ZIPBALL ; \
    rm -f $RPM; \
    mkdir -p $TOOLS/$TOOL_NAME/$TOOL_VERSION/bin ; \
    ln -s /usr/local/bin/bcl2fastq $TOOLS/$TOOL_NAME/$TOOL_VERSION/bin/bcl2fastq ; \
    ln -s $TOOLS/$TOOL_NAME/$TOOL_VERSION $TOOLS/$TOOL_NAME/current ;



############
# BEDTOOLS #
############

ENV TOOL_NAME=bedtools
ENV TOOL_VERSION=2.27.1
ENV TARBALL_LOCATION=https://github.com/arq5x/bedtools2/releases/download/v$TOOL_VERSION
ENV TARBALL=$TOOL_NAME-$TOOL_VERSION.tar.gz
ENV TARBALL_FOLDER=bedtools2
ENV DEST=$TOOLS/$TOOL_NAME/$TOOL_VERSION
ENV PATH=$TOOLS/$TOOL_NAME/$TOOL_VERSION/bin:$PATH

RUN wget $TARBALL_LOCATION/$TARBALL ; \
    tar xf $TARBALL ; \
    rm -rf $TARBALL ; \
    cd $TARBALL_FOLDER ; \
    make prefix=$TOOLS/$TOOL_NAME/$TOOL_VERSION install ; \
    cd ../ ; \
    rm -rf $TARBALL_FOLDER ; \
    ln -s $TOOLS/$TOOL_NAME/$TOOL_VERSION $TOOLS/$TOOL_NAME/current ;



###########
# BOWTIE2 #
###########

ENV TOOL_NAME=bowtie2
ENV TOOL_VERSION=2.3.4.3
ENV GIT=https://github.com/BenLangmead/bowtie2.git
ENV DEST=$TOOLS/$TOOL_NAME/$TOOL_VERSION
ENV PATH=$TOOLS/$TOOL_NAME/$TOOL_VERSION/bin:$PATH

RUN git clone $GIT ; \
    cd $TOOL_NAME ; \
    make -j ; \
    mkdir -p $DEST/bin ; \
    cp LICENSE MANUAL TUTORIAL VERSION NEWS AUTHORS $DEST ; \
    cp bowtie* $DEST/bin ; \
    rm -rf .git ; \
    cd ../ ;  \
    rm -rf $TOOL_NAME ; \
    ln -s $TOOLS/$TOOL_NAME/$TOOL_VERSION/ $TOOLS/$TOOL_NAME/current ;



#######
# BWA #
#######

ENV TOOL_NAME=bwa
ENV TOOL_VERSION=0.7.17
ENV TARBALL_NAME=bwa
ENV TARBALL_LOCATION=https://sourceforge.net/projects/bio-bwa/files/bwa-$TOOL_VERSION.tar.bz2/download
ENV TARBALL=$TARBALL_NAME-$TOOL_VERSION.tar.bz2
ENV DEST=$TOOLS/$TOOL_NAME/$TOOL_VERSION
ENV PATH=$TOOLS/$TOOL_NAME/$TOOL_VERSION/bin:$PATH

RUN wget $TARBALL_LOCATION -O $TARBALL ; \
    tar xf $TARBALL ; \
    rm -rf $TARBALL ; \
    cd $TARBALL_NAME-$TOOL_VERSION ; \
    make ; \
    mkdir -p $TOOLS/$TOOL_NAME/$TOOL_VERSION/bin ; \
    cp bwa $TOOLS/$TOOL_NAME/$TOOL_VERSION/bin ; \
    cd ../ ; \
    rm -rf $TARBALL_NAME-$TOOL_VERSION ; \
    ln -s $TOOLS/$TOOL_NAME/$TOOL_VERSION/ $TOOLS/$TOOL_NAME/current ;



##########
# FASTQC #
##########

ENV TOOL_NAME=fastqc
ENV TOOL_VERSION=0.11.8
ENV TARBALL_NAME=FastQC
ENV TARBALL_LOCATION=http://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v$TOOL_VERSION.zip
ENV TARBALL=fastqc_v$TOOL_VERSION.zip
ENV DEST=$TOOLS/$TOOL_NAME/$TOOL_VERSION
ENV PATH=$TOOLS/$TOOL_NAME/$TOOL_VERSION/bin:$PATH

RUN wget $TARBALL_LOCATION -O $TARBALL ; \
    unzip $TARBALL ; \
    rm -rf $TARBALL ; \
    cd $TARBALL_NAME ; \
    mkdir -p $TOOLS/$TOOL_NAME/$TOOL_VERSION/bin ; \
    chmod u+x fastqc ; \
    cp -R * $TOOLS/$TOOL_NAME/$TOOL_VERSION/bin ; \
    cd ../ ; \
	rm -rf $TARBALL_NAME ; \
    ln -s $TOOLS/$TOOL_NAME/$TOOL_VERSION/ $TOOLS/$TOOL_NAME/current ;



##########
# FATBAM #
##########

ENV TOOL_NAME=fatbam
ENV TOOL_VERSION=0.9.9b
ENV TARBALL_LOCATION=https://gitlab.bioinfo-diag.fr/Strasbourg/FATBAM/repository/0.9.9b
ENV TARBALL=archive.tar.gz
ENV TARBALL_FOLDER=archive
ENV DEST=$TOOLS/$TOOL_NAME/$TOOL_VERSION
ENV PATH=$TOOLS/$TOOL_NAME/$TOOL_VERSION/bin:$PATH

RUN wget $TARBALL_LOCATION/$TARBALL ; \
    tar xf $TARBALL ; \
    rm -rf $TARBALL ; \
    mkdir -p $TOOLS/$TOOL_NAME/$TOOL_VERSION/ ; \
    cp $(ls ${TOOL_NAME^^}-$TOOL_VERSION* -d)/* $TOOLS/$TOOL_NAME/$TOOL_VERSION/ -R ; \
    rm -rf $(ls ${TOOL_NAME^^}-$TOOL_VERSION* -d) ; \
    ln -s $TOOLS/$TOOL_NAME/$TOOL_VERSION $TOOLS/$TOOL_NAME/current ; \
    chmod 0775 $TOOLS/$TOOL_NAME/$TOOL_VERSION $TOOLS/$TOOL_NAME/current -R ; \
    ln -s $TOOLS/$TOOL_NAME/$TOOL_VERSION/ $TOOLS/$TOOL_NAME/current ;



########
# GATK #
########

ENV TOOL_NAME=gatk
ENV TOOL_VERSION=3.8-1-0
ENV TARBALL_RELEASE=gf15c1c3ef
ENV TARBALL_NAME=GenomeAnalysisTK-$TOOL_VERSION-$TARBALL_RELEASE
ENV TARBALL_LOCATION=https://software.broadinstitute.org/gatk/download/auth?package=GATK-archive&version=$TOOL_VERSION-$TARBALL_RELEASE
ENV TARBALL=GenomeAnalysisTK-$TOOL_VERSION.tar.bz2
ENV DEST=$TOOLS/$TOOL_NAME/$TOOL_VERSION
ENV PATH=$TOOLS/$TOOL_NAME/$TOOL_VERSION/bin:$PATH

RUN wget $TARBALL_LOCATION -O $TARBALL ; \
    tar -xf $TARBALL ; \
    rm -rf $TARBALL ; \
    cd $TARBALL_NAME ; \
    mkdir -p $TOOLS/$TOOL_NAME/$TOOL_VERSION/bin ; \
    cp GenomeAnalysisTK.jar $TOOLS/$TOOL_NAME/$TOOL_VERSION/bin ; \
	cd ../ ; \
	rm -rf $TARBALL_NAME ; \
    ln -s $TOOLS/$TOOL_NAME/$TOOL_VERSION/ $TOOLS/$TOOL_NAME/current ;



###########
# HOWARD #
###########

ENV TOOL_NAME=howard
ENV TOOL_VERSION=0.9.14b
ENV TARBALL_LOCATION=https://gitlab.bioinfo-diag.fr/Strasbourg/HOWARD/repository/$TOOL_VERSION
ENV TARBALL=archive.tar.gz
ENV TARBALL_FOLDER=archive
ENV TOOL_DATABASE_FOLDER=/home/TOOLS/databases
ENV DEST=$TOOLS/$TOOL_NAME/$TOOL_VERSION
ENV PATH=$TOOLS/$TOOL_NAME/$TOOL_VERSION/bin:$PATH

RUN wget $TARBALL_LOCATION/$TARBALL ; \
    tar xf $TARBALL ; \
    rm -rf $TARBALL ; \
    mkdir -p $TOOLS/$TOOL_NAME/$TOOL_VERSION/ ; \
    cp $(ls ${TOOL_NAME^^}-$TOOL_VERSION* -d)/* $TOOLS/$TOOL_NAME/$TOOL_VERSION/ -R ; \
    rm -rf $(ls ${TOOL_NAME^^}-$TOOL_VERSION* -d) ; \
    ln -s $TOOLS/$TOOL_NAME/$TOOL_VERSION $TOOLS/$TOOL_NAME/current ; \
    chmod 0775 $TOOLS/$TOOL_NAME/$TOOL_VERSION $TOOLS/$TOOL_NAME/current -R ; \
	mkdir -p $DATABASES ; \
	ln -s $DATABASES $TOOL_DATABASE_FOLDER ; \
	ln -s $TOOLS/$TOOL_NAME/$TOOL_VERSION $TOOLS/$TOOL_NAME/current ;



############
# IGVTOOLS #
############

ENV TOOL_NAME=igvtools
ENV TOOL_MAIN_VERSION=2.4
ENV TOOL_VERSION=$TOOL_MAIN_VERSION.16
ENV TARBALL_NAME=igvtools_$TOOL_VERSION
ENV TARBALL=$TARBALL_NAME.zip
ENV TARBALL_FOLDER=IGVTools
ENV TARBALL_LOCATION=http://data.broadinstitute.org/igv/projects/downloads/$TOOL_MAIN_VERSION/$TARBALL
ENV DEST=$TOOLS/$TOOL_NAME/$TOOL_VERSION
ENV PATH=$TOOLS/$TOOL_NAME/$TOOL_VERSION/bin:$PATH

RUN wget $TARBALL_LOCATION -O $TARBALL ; \
    unzip $TARBALL ; \
    rm -rf $TARBALL ; \
    cd $TARBALL_FOLDER ; \
    mkdir -p $TOOLS/$TOOL_NAME/$TOOL_VERSION/bin ; \
    cp -R * $TOOLS/$TOOL_NAME/$TOOL_VERSION/bin ; \
	cd ../ ; \
	rm -rf $TARBALL_FOLDER ; \
    ln -s $TOOLS/$TOOL_NAME/$TOOL_VERSION/ $TOOLS/$TOOL_NAME/current ;



##########
# ITDSEEK #
###########

ENV TOOL_NAME=itdseek
ENV TOOL_VERSION=1.2-2
ENV TARBALL_LOCATION=https://github.com/tommyau/itdseek/zipball/master
ENV TARBALL=$TOOL_NAME-$TOOL_VERSION.zip
ENV TARBALL_FOLDER=tommyau-itdseek-*
ENV DEST=$TOOLS/$TOOL_NAME/$TOOL_VERSION
ENV PATH=$TOOLS/$TOOL_NAME/$TOOL_VERSION/bin:$PATH

RUN wget $TARBALL_LOCATION/$TARBALL -O $TARBALL ; \
    unzip $TARBALL ; \
    rm -rf $TARBALL ; \
	mkdir -p $TOOLS/$TOOL_NAME/$TOOL_VERSION/bin ; \
	cp $TARBALL_FOLDER/* $TOOLS/$TOOL_NAME/$TOOL_VERSION/bin ; \
    rm -rf $TARBALL_FOLDER ; \
    ln -s $TOOLS/$TOOL_NAME/$TOOL_VERSION/ $TOOLS/$TOOL_NAME/current ;



########
# JAVA #
########

ENV TOOL_NAME=java
ENV TOOL_VERSION=current
RUN mkdir -p $TOOLS/$TOOL_NAME/$TOOL_VERSION/bin && \
	ln -s /usr/bin/java $TOOLS/$TOOL_NAME/$TOOL_VERSION/bin/java ;



##########
# PICARD #
##########

ENV TOOL_NAME=picard
ENV TOOL_VERSION=2.18.5
ENV JAR_LOCATION=https://github.com/broadinstitute/picard/releases/download/$TOOL_VERSION
ENV JAR=picard.jar
ENV DEST=$TOOLS/$TOOL_NAME/$TOOL_VERSION
ENV PATH="$TOOLS/$TOOL_NAME/$TOOL_VERSION/bin:$PATH"

RUN wget $JAR_LOCATION/$JAR ; \
    mkdir -p $TOOLS/$TOOL_NAME/$TOOL_VERSION/bin ; \
    mv $JAR $TOOLS/$TOOL_NAME/$TOOL_VERSION/bin ; \
    ln -s $TOOLS/$TOOL_NAME/$TOOL_VERSION $TOOLS/$TOOL_NAME/current ;



############
# SAMTOOLS #
############

ENV TOOL_NAME=samtools
ENV TOOL_VERSION=1.8
ENV TARBALL_LOCATION=https://github.com/samtools/$TOOL_NAME/releases/download/$TOOL_VERSION/
ENV TARBALL=$TOOL_NAME-$TOOL_VERSION.tar.bz2
ENV DEST=$TOOLS/$TOOL_NAME/$TOOL_VERSION
ENV PATH=$TOOLS/$TOOL_NAME/$TOOL_VERSION/bin:$PATH

RUN wget $TARBALL_LOCATION/$TARBALL ; \
    tar xf $TARBALL ; \
    rm -rf $TARBALL ; \
    cd $TOOL_NAME-$TOOL_VERSION ; \
    make prefix=$TOOLS/$TOOL_NAME/$TOOL_VERSION install ; \
    cd ../ ; \
    rm -rf $TOOL_NAME-$TOOL_VERSION ; \
    ln -s $TOOLS/$TOOL_NAME/$TOOL_VERSION $TOOLS/$TOOL_NAME/current ;



##########
# SNPEFF #
##########

ENV TOOL_NAME=snpeff
ENV TOOL_VERSION=4.3t
ENV TOOL_VERSION_FOR_FILE=4_3t
ENV TARBALL_LOCATION=https://sourceforge.net/projects/snpeff/files
ENV TARBALL="snpEff_v"$TOOL_VERSION_FOR_FILE"_core.zip"
ENV TARBALL_FOLDER=snpeff_folder
ENV TOOL_DATABASE_FOLDER=$TOOLS/$TOOL_NAME/$TOOL_VERSION/bin/data
ENV TOOL_DATABASE_FOLDER_LINK=$DATABASES/snpeff_sources/4.3t
ENV DEST=$TOOLS/$TOOL_NAME/$TOOL_VERSION
ENV PATH=$TOOLS/$TOOL_NAME/$TOOL_VERSION/bin:$PATH

RUN wget $TARBALL_LOCATION/$TARBALL ; \
    unzip $TARBALL -d $TARBALL_FOLDER ; \
	rm $TARBALL ; \
    mkdir -p $TOOLS/$TOOL_NAME/$TOOL_VERSION/bin/ ; \
    cp $TARBALL_FOLDER/*/*jar $TOOLS/$TOOL_NAME/$TOOL_VERSION/bin/ -R ; \
    cp $TARBALL_FOLDER/*/*.config $TOOLS/$TOOL_NAME/$TOOL_VERSION/bin/ -R ; \
    rm -rf $TARBALL_FOLDER ; \
	ln -s $TOOLS/$TOOL_NAME/$TOOL_VERSION $TOOLS/$TOOL_NAME/current ; \
	mkdir -p $TOOL_DATABASE_FOLDER_LINK ; \
	mkdir -p $TOOL_DATABASE_FOLDER ; \
	ln -s $TOOL_DATABASE_FOLDER_LINK $TOOL_DATABASE_FOLDER ;


###########
# VARSCAN #
###########

ENV TOOL_NAME=varscan
ENV TOOL_VERSION=2.4.3
ENV TARBALL_RELEASE=
ENV TARBALL_NAME=VarScan.v$TOOL_VERSION.jar
ENV TARBALL_LOCATION=https://github.com/dkoboldt/varscan/raw/master/$TARBALL_NAME
ENV TARBALL=VarScan.jar
ENV DEST=$TOOLS/$TOOL_NAME/$TOOL_VERSION
ENV PATH=$TOOLS/$TOOL_NAME/$TOOL_VERSION/bin:$PATH

RUN wget $TARBALL_LOCATION -O $TARBALL ; \
    mkdir -p $TOOLS/$TOOL_NAME/$TOOL_VERSION/bin ; \
	cp $TARBALL $TOOLS/$TOOL_NAME/$TOOL_VERSION/bin ; \
	rm $TARBALL ; \
    ln -s $TOOLS/$TOOL_NAME/$TOOL_VERSION/ $TOOLS/$TOOL_NAME/current ;



############
# VCFTOOLS #
############

ENV TOOL_NAME=vcftools
ENV TOOL_VERSION=0.1.14
ENV TARBALL_NAME=vcftools
ENV TARBALL_LOCATION=https://github.com/samtools/$TARBALL_NAME/releases/download/$TOOL_VERSION/
ENV TARBALL=$TARBALL_NAME-$TOOL_VERSION.tar.bz2
ENV DEST=$TOOLS/$TOOL_NAME/$TOOL_VERSION
ENV PATH=$TOOLS/$TOOL_NAME/$TOOL_VERSION/bin:$PATH

RUN git clone https://github.com/cgrlab/vcftools.git ; \
    cd vcftools; \
    git checkout tags/v0.1.14 ; \
    ./autogen.sh ; \
    ./configure ; \
    make prefix=$TOOLS/$TOOL_NAME/$TOOL_VERSION ; \
    make prefix=$TOOLS/$TOOL_NAME/$TOOL_VERSION install ; \
    cd ../ ; \
	rm -rf TARBALL_NAME ; \
    ln -s $TOOLS/$TOOL_NAME/$TOOL_VERSION/ $TOOLS/$TOOL_NAME/current ;



#########
# LATEX #
#########

ENV TOOL_NAME=texlive
ENV TOOL_VERSION=current
ENV TARBALL_LOCATION=http://mirror.ctan.org/systems/texlive/tlnet
ENV TARBALL=install-tl-unx.tar.gz
ENV TARBALL_FOLDER=install-tl-*
ENV DEST=$TOOLS/$TOOL_NAME/$TOOL_VERSION
ENV PATH=$TOOLS/$TOOL_NAME/$TOOL_VERSION/bin:$PATH


RUN wget $TARBALL_LOCATION/$TARBALL ; \
    tar xf $TARBALL ; \
    rm -rf $TARBALL ; \
    cd $TARBALL_FOLDER ; \
    mkdir -p $TOOLS/$TOOL_NAME/ ; \
    echo 'I' | ./install-tl -scheme=scheme-basic ; \
    cd ../ ; \
    rm -rf $TARBALL_FOLDER ; \
	mkdir -p $TOOLS/$TOOL_NAME/$(basename $(ls -d /usr/local/texlive/2*)) ; \
    ln -s /usr/local/texlive/$(basename $(ls -d /usr/local/texlive/2*))/bin/x86_64-linux/ $TOOLS/$TOOL_NAME/$(basename $(ls -d /usr/local/texlive/2*))/bin ; \
	ln -s $TOOLS/$TOOL_NAME/$(basename $(ls -d /usr/local/texlive/2*))/ $TOOLS/$TOOL_NAME/current ; \
	$TOOLS/$TOOL_NAME/current/bin/tlmgr install tcolorbox ; \
	$TOOLS/$TOOL_NAME/current/bin/tlmgr install trimspaces ; \
	$TOOLS/$TOOL_NAME/current/bin/tlmgr install lastpage ; \
	$TOOLS/$TOOL_NAME/current/bin/tlmgr install numprint ; \
	$TOOLS/$TOOL_NAME/current/bin/tlmgr install multirow ; \
	$TOOLS/$TOOL_NAME/current/bin/tlmgr install supertabular ; \
	$TOOLS/$TOOL_NAME/current/bin/tlmgr install lipsum ; \
	$TOOLS/$TOOL_NAME/current/bin/tlmgr install siunitx ; \
	$TOOLS/$TOOL_NAME/current/bin/tlmgr install tocloft ; \
	$TOOLS/$TOOL_NAME/current/bin/tlmgr install l3kernel ; \
	$TOOLS/$TOOL_NAME/current/bin/tlmgr install l3packages ; \
	$TOOLS/$TOOL_NAME/current/bin/tlmgr install fancybox ; \
	$TOOLS/$TOOL_NAME/current/bin/tlmgr install listings ; \
	$TOOLS/$TOOL_NAME/current/bin/tlmgr install xcolor ; \
	$TOOLS/$TOOL_NAME/current/bin/tlmgr install pgf ; \
	$TOOLS/$TOOL_NAME/current/bin/tlmgr install collcell ; \
	$TOOLS/$TOOL_NAME/current/bin/tlmgr install etoolbox ; \
	$TOOLS/$TOOL_NAME/current/bin/tlmgr install environ ; \
	$TOOLS/$TOOL_NAME/current/bin/tlmgr install caption ; \
	$TOOLS/$TOOL_NAME/current/bin/tlmgr install subfig ; \
	$TOOLS/$TOOL_NAME/current/bin/tlmgr install float ; \
	$TOOLS/$TOOL_NAME/current/bin/tlmgr install setspace ; \
	$TOOLS/$TOOL_NAME/current/bin/tlmgr install fp ;





######################
# YUM REMOVE & CLEAR #
######################

RUN yum erase -y $YUM_REMOVE ; yum clean all ;



##############################
# WORKDIR / ENTRYPOINT / CMD #
##############################


WORKDIR "/STARK/data"

#ENTRYPOINT [ "/bin/bash" ]
# find /Users/lebechea/Documents/NGS/Dockerfile -maxdepth 2 -mindepth 2 -type d
#CMD [ "ls", "-l", "/home/TOOLS/tools" ]
CMD [ "find", "/STARK/tools", "-maxdepth", "2", "-mindepth", "2", "-type", "d" ]
