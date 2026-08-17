"""Python FDICAパッケージの公開窓口。

利用者が ``from python_fdica import bssAuxFdica`` と簡潔に読み込めるよう、
内部モジュールの主関数をパッケージ直下へ公開する。
"""

from .bssAuxFdica import bssAuxFdica

# ワイルドカードimportを行った場合に公開する名前を明示する。
__all__ = ["bssAuxFdica"]
