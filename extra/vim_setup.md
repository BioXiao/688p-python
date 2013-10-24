A Simple Vim Setup
==================
Below are some basic steps you can take to improve the appearance and
functionality of Vim.

First, install the [Pathogen](https://github.com/tpope/vim-pathogen) plugin 
manager for Vim. There are numerous plugins available for Vim to extend its
basic functionality. Pathogen helps to make the process of adding and removing
those plugins simpler.

```
mkdir -p ~/.vim/autoload ~/.vim/bundle; \
curl -Sso ~/.vim/autoload/pathogen.vim \
https://raw.github.com/tpope/vim-pathogen/master/autoload/pathogen.vim
```

Next, install [sensible.vim](https://github.com/tpope/vim-sensible). 
sensible.vim is a small plugin written by the same author of Pathogen which
simply applied some "sensible" default settings to Vim. This will fix things
like issues with the backspace key, indentation, etc.

```
cd ~/.vim/bundle
git clone git://github.com/tpope/vim-sensible.git
```

Install [vim-airline](https://github.com/bling/vim-airline) for a more obvious
indication of when you are working in INSERT/NORMAL mode.

```
git clone https://github.com/bling/vim-airline
```

There are many themes available for Vim. One of the ones I like is called
[Molokai](https://github.com/tomasr/molokai):

```
git clone https://github.com/tomasr/molokai.git
```

In order for these changes to take effect, let's make a few changes to
your vim configuration file ([.vimrc](http://vimdoc.sourceforge.net/htmldoc/starting.html))
Using vim, editing ~/.vimrc and add the following lines:

```
" Load pathogen
execute pathogen#infect()

" Enable syntax highlighting
set background=dark                                                                                     
colorscheme molokai
syntax on

" Other
set cursorline
set number
```

Finally, if you are using Grace at UMD, you will want to tell your shell to
use /usr/bin/vim (Vim 7.2) instead of /usr/local/bin/vim (Vim 7.0). To do this,
run the following command:

```
echo 'alias vim=/usr/bin/vim' >> .bashrc.mine
```
Done!
