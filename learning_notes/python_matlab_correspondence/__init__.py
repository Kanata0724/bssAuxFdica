# 対応注: 学習用コピー：元ファイル python_fdica/__init__.py（learning / b934fef 時点）。
# 対応注: パッケージ公開設定はPython独自。MATLABに同じ構文の対応箇所はない。
"""Python FDICAパッケージの公開窓口。

利用者が ``from python_fdica import bssAuxFdica`` と簡潔に読み込めるよう、
内部モジュールの主関数をパッケージ直下へ公開する。
"""

from .bssAuxFdica import bssAuxFdica

# ワイルドカードimportを行った場合に公開する名前を明示する。
__all__ = ["bssAuxFdica"]
