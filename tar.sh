#!/bin/bash 
cd $HOME/code/eampa/
tar cvzf eampa.tar.gz * --directory="$HOME/code/eampa" --exclude="$HOME/code/eampa/.git" \
 --exclude="eampa.tar.gz"