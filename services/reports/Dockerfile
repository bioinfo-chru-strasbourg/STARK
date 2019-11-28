FROM ubuntu:16.04

# sudo docker run -d -p 8080:80 -v '/path/to/logs':/var/log/apache2 -v '/path/to/site/files':/var/www/html image_name
# docker run -p 4299:80 -v /Users/lebechea/Documents/NGS/GIT/Strasbourg/STARK/service/reports:/var/www/html -v /Users/lebechea/Documents/NGS/DATA/REPOSITORY:/var/www/html/repositories/REPOSITORY --name=STARK_REPORTS stark_httpd
# docker run -p 4298:80 -v /Users/lebechea/Documents/NGS/GIT/Strasbourg/STARK/service/reports:/var/www/html -v /Users/lebechea/Documents/NGS/DATA/REPOSITORY:/var/www/html/repositories/REPOSITORY --entrypoint=bash -ti stark_httpd

EXPOSE 80

RUN apt-get update
RUN apt-get -y install apache2
RUN apt-get -y install python
RUN apt-get -y install python-pil python-matplotlib python-scipy python-sklearn
RUN apt-get -y install php libapache2-mod-php

#load apache cgi module
RUN a2enmod cgi
RUN service apache2 restart

COPY index.php /var/www/html/index.php
COPY assets /var/www/html/assets

#enable cgi in the website root
#second block to allow .htaccess
RUN echo "                          \n \
<Directory /var/www/html>           \n \
   Options +ExecCGI                 \n \
   AddHandler cgi-script .php .py   \n \
   DirectoryIndex index.php         \n \
</Directory>                        \n \
" >> /etc/apache2/apache2.conf

RUN chmod -R u+rwx,g+x,o+x /var/www/html

RUN ln -sf /usr/bin/python /usr/local/bin/python

CMD /usr/sbin/apache2ctl -D FOREGROUND
