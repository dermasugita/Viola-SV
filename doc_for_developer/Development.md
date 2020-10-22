# PySgt(Structural Genome Translator on Python)開発用ドキュメント

## 推奨開発環境

### OSに関して
client OS: Any<br>
remote OS: Linux

### Pythonのバージョン
- pyenvでの管理がお勧め

python3.6, 3.7, 3.8で動作保証したいです。
テストを3つのバージョンで行う必要があり、各バージョンのインタープリタへのPathを同時に通しておく必要があるためです。

インストールは[こちら](https://github.com/pyenv/pyenv)

プロジェクトディレクトリで<br>`pyenv local 3.6.? 3.7.? 3.8.?`<br>とすることで、全てのバージョンにPathが通る。?に自分がインストールしたマイナーバージョンを指定。

### パッケージ管理
- **pipenv**または**virtualenv**でのパッケージ管理を推奨

パッケージ管理はテスト環境構築を目的とします。
pyenvに直接テストパッケージをインストールすることはお勧めできないです。
自分はpipenvを利用しています。管理が楽チンでおすすめです。

pipenvを使う場合、以下のように設定します(インストールしているのはテストに必要なものです)。
```shell
pip instal pipenv
pipenv --python 3
pipenv install pytest tox
```
pipenvを有効にしたい場合は、`pipenv shell`を入力。

詳しい説明は[この辺](https://qiita.com/y-tsutsu/items/54c10e0b2c6b565c887a)とかみてください。


## テスト環境
テストは
- pytest
- tox

で行います。

### pytestに関して
**pytest**はPythonのテストツールです。pythonのバージョンにこだわらなければ、自分のpipenv環境に開発モードで自分のパッケージをインストールし、**pytest**だけでテスト環境を構築することもできます。
```
$ pytest [オプションとか引数とか]
```
というコマンドでテストコードを実行できます。詳しくは[ここ](https://qiita.com/everylittle/items/1a2748e443d8282c94b2)とかで。

### toxに関して
今回は3つのpythonバージョンでテストしたいため、**tox**というパッケージを利用します。

**tox**はpythonのバージョンごとに、`.tox/`以下に自動でvirtualenvを作成し、それぞれテストコードが正常に動作するかをテストしてくれます。`.tox/`以下のvirtualenvは、パッケージの`setup.py`を元に依存パッケージを自動でインストールしてくれます。

`/path/to/the/project_root/tox.ini`に以下の設定をします。

`file tox.ini`の内容↓
```
[tox]
envlist = py36,py37,py38

[testenv]
deps = pytest
commands = pytest {posargs}
```

これで`tox`コマンドだけでテストができます。pytestにオプションや引数を渡したければ、

```
tox -- 何かオプション
```
でできます。詳しくは[こちら](https://tox.readthedocs.io/en/latest/example/pytest.html)。


## コーディングの決まり
以下のコーディング規則でお願いします。途中途中で決めているので、僕自身の初期のコードのなかにも遵守されていないものがあるかもしれません。気がついたら直しておいてください（怖い場合は連絡ください）。
### 命名規則
- 変数のprefix
<br>`list dictionary set tuple ndarray pd.Series pd.DataFrame`の場合は、それぞれ`ls_ dict_ set_ tuple_ arr_ ser_ df_`というprefixをつけてください。 intとかstrとかはテキトーでOKです。
- 変数、関数、メソッドはアンダースコア区切りで全て小文字
- クラス、例外、自定義型は最初大文字+大文字区切り
- 定数はアンダースコア区切りで全て大文字
- クラス内変数・メソッドはアンダースコアで始めてください。 例: `_inner_class_method()`

### Docstringsを書きましょう
- [numpy style](https://numpydoc.readthedocs.io/en/latest/format.html)でお願いします。
- classのDocstringsはclass定義の真下に書き、`__init__()`の下には何も書かないこと。

### 型ヒントは出来るだけ書いといてください
- 変数の引数や出力の型が定まっている場合は、型ヒントを書いておいてください
- 確定できない変数は、悩んで時間をとられるよりは何も書かずに次に進んで良いです
- 新しい型を定義したい場合は、`_typing.py`に追加してください。

### Errorの定義
- Exceptionを自分で定義したい場合は`_exception.py`に書いてください。

## Gitでの管理に関して
### 基本事項
- コードの管理はGitHubで行います。
- 開発のメインはdevelopブランチです。masterには触らないこと。
- developは直接編集しないこと。
- developからfeature/somethingという名前のブランチを切って作業してください。

### 最新のorigin/developを自分の作業ブランチ(feature/something)に取り込みたい
- `git pull --rebase origin develop`を叩いてください。
- `git push origin feature/something`が通らない場合は、`--force-with-lease`をつけてください。

### developへのマージ
- developへのマージはGitHub上で行います
- `git pull --rebase origin develop`を叩いてから`git push --force-with-lease origin feature/something`し、GitHub上でPull Requestしてください

## こうしておけば基本的に間違いない、という流れ

### 自分用のブランチを切る
- まず`git checkout develop`でdevelopに行く
- `git checkout -b feature/something`

### 自分のブランチでの作業してリモートにpushするまでの流れ
- `git checkout feature/something`で自分のブランチへ行く
- 今いるブランチは`git branch`でわかる
- コーディングする
- 一通りコーディングできたら`git add .`してコミット`git commit -m "some message"`
- pushする際は、まずdevelopをrebaseする。`git pull --rebase origin develop`
    - conflictした際は手動で解消し、addする。`git add コンフリクト解消したファイル`
    - `git rebase --continue`
    - 下に合流
- rebaseできたらリモートにプッシュ(初めて) `git push origin feature/something` 
- プッシュ2回目以降 `git push --force-with-lease origin feature/something`

### 結構ちゃんとできたのでdevelopにマージしたいときの流れ
- 真上のフローで`git push --force-with-lease origin feature/something`までやる
- GitHub上でPull Requestする
- 皆にSlackする
- 杉田が確認しマージする

### やってしまった！
- まずは`git log`して落ち着く(多くの場合commit idから復活できる)
- 対応がわからなければslackで


## その他
- 規則を破った人を怒るのはやめましょう。
- 困ったらslackで
- ディレクトリ構成はpandasを参考にしています。