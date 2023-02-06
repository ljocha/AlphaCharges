FROM ubuntu:22.04

RUN apt update -y && apt upgrade -y
RUN apt-get --yes install python3-pip python3-venv apache2 libapache2-mod-wsgi-py3 gemmi
RUN apt install -y vim

WORKDIR /opt
RUN python3 -m venv venv

ARG activate=". /opt/venv/bin/activate"
RUN ${activate} && pip install flask requests scipy numba numpy pdb2pqr scikit-learn rdkit==2022.03.5

RUN mkdir -p /opt/AlphaCharges/app
COPY app/ /opt/AlphaCharges/app

# RUN chown -R ubuntu:ubuntu /opt
RUN rm -f /etc/apache2/sites-available/*
COPY AlphaCharges.conf /etc/apache2/sites-available/
RUN chown -R www-data:www-data /opt
RUN chmod o+rx AlphaCharges/app/AlphaCharges.wsgi
RUN chmod o+rx AlphaCharges/app/routes.py
RUN a2ensite AlphaCharges.conf
RUN a2enmod ssl
RUN a2enmod brotli
RUN a2enmod http2


# sudo systemctl restart apache2
EXPOSE 80
EXPOSE 443

ENTRYPOINT apachectl -D FOREGROUND
