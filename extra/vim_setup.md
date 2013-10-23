A Simple Vim Setup
==================
Below are some basic steps you can take to improve the appearance and
functionality of Vim.

1. Install [Pathogen](https://github.com/tpope/vim-pathogen) plugin manager for
vim

```
mkdir -p ~/.vim/autoload ~/.vim/bundle; \
curl -Sso ~/.vim/autoload/pathogen.vim \
https://raw.github.com/tpope/vim-pathogen/master/autoload/pathogen.vim
```

2. Install [sensible.vim](https://github.com/tpope/vim-sensible)

```
cd ~/.vim/bundle
git clone git://github.com/tpope/vim-sensible.git
```

3. Install [vim-airline](https://github.com/bling/vim-airline)

```
git clone https://github.com/bling/vim-airline ~/.vim/bundle/vim-airline
```

4. Install Molokai theme

```
mkdir ~/.vim/colors
cd ~/.vim/colors
wget https://github.com/tomasr/molokai/raw/master/colors/molokai.vim
```

5. Install solarized theme

```
cd ~/.vim/bundle
git clone git://github.com/altercation/vim-colors-solarized.git
```

4. Edit your [.vimrc](http://vimdoc.sourceforge.net/htmldoc/starting.html) to 
include:

```
" Load pathogen
execute pathogen#infect()

" Enable syntax highlighting
set background=dark                                                                                     
colorscheme molokai
syntax on

" Other
set cursorline
```
