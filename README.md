# Hirschberg's algorithm
Hirschberg's algorithm is a dynamic programming algorithm that finds the optimal sequence alignment between two strings **X** with length *n* and **Y** with length *m*.
In comparison with Needleman-Wunsch algorithm, Hirschberg's algorithm is more space efficient:

| Algorithm | Time complexity | Space complexity |
|---|---|---|
| Needleman-Wunsch | *O(n \* m)* | *O(n \* m)* |
| Hirschberg | *O(n \* m)* | *O(n + m)* |

### Task
* Implement Hirschberg's algorithm - you can use a template `hirschberg_template.R`


More information can be also file on [wikipedia](https://en.wikipedia.org/wiki/Hirschberg%27s_algorithm).

<details>
<summary>Spoilers! Help</summary>

#### Pseudo code of Hirschberg's algorithm
```
function Hirschberg(X, Y)
    Z = ""
    W = ""
    if length(X) == 0
        for i=1 to length(Y)
            Z = Z + '-'
            W = W + Y(i)
        end
    else if length(Y) == 0
        for i=1 to length(X)
            Z = Z + X(i)
            W = W + '-'
        end
    else if length(X) == 1 or length(Y) == 1
        (Z, W) = NeedlemanWunsch(X, Y)
    else
        xlen = length(X)
        xmid = length(X) / 2
        ylen = length(Y)

        ScoreL = NWScore(X(1:xmid), Y)
        ScoreR = NWScore(rev(X(xmid+1:xlen)), rev(Y))
        ymid = arg max ScoreL + rev(ScoreR)

        (Z, W) = Hirschberg(X(1:xmid), Y(1:ymid)) + Hirschberg(X(xmid+1:xlen), Y(ymid+1:ylen))
    end
    return (Z, W)
```
</details>
<details>
<summary>Download files from GitHub</summary>
<details>
<summary>Git settings</summary>

> * Configure the Git editor
> ```bash
> git config --global core.editor notepad
> ```
> * Configure your name and email address
> ```bash
> git config --global user.name "Zuzana Nova"
> git config --global user.email z.nova@vut.cz
> ```
> * Check current settings
> ```bash
> git config --global --list
> ```
>
</details>

* Create a fork on your GitHub account. 
  On the GitHub page of this repository find a <kbd>Fork</kbd> button in the upper right corner.
  
* Cloned forked repository from your GitHub page to a folder in your computer:
```bash
git clone <fork repository address>
```
* In a local repository, set new remote for project repository:
```bash
git remote add upstream https://github.com/mpa-prg/exercise_04.git
```

### Send files to GitHub
Create a new commit and send new changes to your remote repository.
* Add file to a new commit.
```bash
git add <file_name>
```
* Create commit, enter commit message, save the file and close it.
```bash
git commit
```
* Send new commit to your GitHub repository.
```bash
git push origin master
```
</details>
