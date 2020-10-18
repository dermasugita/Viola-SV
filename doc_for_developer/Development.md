# PySgt(Structural Genome Translator on Python)開発用ドキュメント

## 推奨開発環境

### OSに関して
client OS: Any
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

### 変数の命名規則
以下のコーディング規則でお願いします。途中途中で決めているので、僕自身の初期のコードのなかにも遵守されていないものがあるかもしれません。気がついたら直しておいてください（怖い場合は連絡ください）。
- 変数のprefix
<br>`list dictionary set tuple ndarray pd.Series pd.DataFrame`の場合は、それぞれ`ls_ dict_ set_ tuple_ arr_ ser_ df_`というprefixをつけてください。 

- 変数、関数、メソッドはアンダースコア区切りで全て小文字
- クラス、例外、自定義型は最初大文字+大文字区切り
- 定数はアンダースコア区切りで全て大文字