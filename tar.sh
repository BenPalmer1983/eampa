#!/bin/bash 
cd $HOME/code/eampa/
tar cvzf eampa.tar.gz * --directory="/code/eampa" \
--exclude="*/.git" \
--exclude="make.sh" \
--exclude="*.docx" \
--exclude="eampa.tar.gz" \
--exclude="tar.sh" \
--exclude="git.sh" \
--exclude="untar.sh" \
--exclude="bin/*.x"