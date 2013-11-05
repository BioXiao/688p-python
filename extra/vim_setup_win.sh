#!/bin/sh
# CBMG 688P: Vim setup script
# Keith Hughitt (khughitt@umd.edu)
# 2013/10/24

# Pathogen
echo "Installing pathogen.vim"
mkdir -p ~/vimfiles/autoload ~/vimfiles/bundle; \
curl -Sso ~/vimfiles/autoload/pathogen.vim \
https://raw.github.com/tpope/vim-pathogen/master/autoload/pathogen.vim

echo "Installing sensible.vim"
cd ~/vimfiles/bundle
git clone git://github.com/tpope/vim-sensible.git

echo "Installing vim-airline"
git clone https://github.com/bling/vim-airline

echo "Installing molokai theme for vim"
git clone https://github.com/tomasr/molokai.git

echo "Creating Vim configuration"
read -d '' vimconf <<-"_EOF_"
" Load pathogen
execute pathogen#infect()

" Enable syntax highlighting
set background=dark
let g:airline_theme = "powerlineish"
colorscheme molokai
syntax on

" Other
set cursorline
set number
_EOF_

echo "$vimconf" >> ~/_vimrc

