# Povu TUI Documentation

## Introduction

Povu TUI is an interactive terminal user interface that helps you explore and analyse variant call data and graphs efficiently without leaving your terminal.

To enter the TUI, run `povu view`.

## Modes

Povu TUI provides two main viewers, each designed for specific tasks:

1. **Graph Explorer**: This viewer allows you to visualize and traverse genomic graphs interactively, making it easier to analyze structural relationships.
2. **VCF Explorer**: This viewer enables browsing, searching, and analyzing variant call format (VCF) data in a tabular format.

## Modes and Command Mode

Press `<SPACE>` to enter command mode.

To switch between the two viewers:
 - From VCF Explorer to Graph Explorer: `<SPACE> + g`
 - From Graph Explorer to VCF Explorer: `<SPACE> + v`

**Navigation**

 - Switch Panes (VCF Explorer): Press `<TAB>` to switch between panes within the VCF Explorer.\
 - Scroll: Use the arrow keys (↑ , ↓, ←, →) to move through content.

**Search**

 - Press `/` to enter search mode. \
In the graph viewer, this enables graph exploration within a 250×250 vertex window. \
To scroll up and down search results use `n` to go forward and `N` to go back.

**Jump to Line**

 - Jump to a specific line: Press `:` and enter the line number.


## VCF Explorer

* Press `<TAB>` to switch between panes within the VCF Explorer.