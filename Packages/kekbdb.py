import oracledb
import rdbtool

_default_connection='readonly/readonly@kekb'

kekbdb=rdbtool.rdb(oracledb.oracledb(_default_connection))
