# PySgt(Structural Genome Translator on Python)開発用ドキュメント

## 推奨開発環境

### OSに関して
client OS: Any
remote OS: Linux

### Pythonのバージョン
python3.6, 3.7, 3.8で動作保証したいです。
テストを3つのバージョンで行う必要があり、各バージョンのインタープリタへのPathを同時に通しておく必要があります。

そのため**pyenv**で管理するのを推奨。

インストールは[こちら](https://github.com/pyenv/pyenv)

プロジェクトディレクトリで<br>`pyenv local 3.6.? 3.7.? 3.8.?`<br>とすることで、全てのバージョンにPathが通る。?に自分がインストールしたマイナーバージョンを指定。

### パッケージ管理
パッケージ管理はテスト環境構築を目的とします。

pyenvに直接テストパッケージをインストールすることはお勧めできないため、virtualenvやpipenvなどが推奨です。
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