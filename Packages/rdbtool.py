import string
import UserList

class _rdbtabcmp :
    def __init__(self,columns) :
        self.columns=columns
        return
    def cmp(self,x,y) :
        for i in self.columns :
            c=cmp(x[i],y[i])
            if c : return c
        return 0

class rdbtab(UserList.UserList) :
    def __init__(self,list,desc=None,sql=None) :
        UserList.UserList.__init__(self,list)
        self.description=desc
        self.sql=sql
        return
    def __setitem__(self,i,item) :
        self.data[i]=item
        self.description=self.sql=None
    def __delitem__(self,i) :
        del self.data[i]
        self.sql=None
    def __getslice__(self,i,j) :
        return self.__class__(self.data[i:j], self.description)
    def __setslice__(self,i,j,list) :
        UserList.UserList.__setslice__(self,i,j,list)
        self.description=self.sql=None
    def __delslice__(self,i,j) :
        del self.data[i:j]
        self.sql=None
    def append(self,item) :
        self.data.append(item)
        self.description=self.sql=None
    def insert(self,i,item) :
        self.data.insert(i,item)
        self.description=self.sql=None
    def remove(self,item) :
        self.data.remove(item)
        self.sql=None
    def reverse(self) :
        self.data.reverse()
        self.sql=None
    def sort(self,*args) :
        if not args :
            self.data.sort()
        elif callable(args[0]) :
            apply(self.data.sort,args)
        else :
            if type(args[0]) not in (type(0),type("")) : args=args[0]
            columns=[]
            for item in args :
               if type(item)==type(0) : columns.append(item)
               else                   : columns.append(self.colindex(item))
            c=_rdbtabcmp(columns)
            self.data.sort(c.cmp)
        self.sql=None
        return
    def colindex(self,name) :
        name=string.upper(name)
        i=0
        for col in self.description :
            if col[0]==name : return i
            i=i+1
        raise ValueError,name+': invalid column name'

class rdb :
    def __init__(self,con=None) :
        self.lastsql=None
        return self.open(con)
    def __del__(self) :
        self.close()
        return
    def __getattr__(self,attr) :
        if not self.con or attr[0]=='_' : raise AttributeError,attr
        try :                   return getattr(self.con,attr)
        except AttributeError : pass
        try :                          return getattr(self.cur,attr)
        except AttributeError,detail : raise AttributeError,detail
    def open(self,con=None) :
        self.con=con
        if con : self.cur=con.cursor()
        return
    def close(self) :
        if not self.con : return
        self.cur.close()
        self.con.close()
        self.con=None
        return
# additional functions
    def select(self,str,params=None) :
        if string.upper(string.split(str,None,1)[0])!='SELECT' :
            str='SELECT '+str
        if params : self.cur.execute(str,params)
        else      : self.cur.execute(str)
        self.lastsql=str
        return self.cur.fetchall()
    def gettab(self,table) :
        return self.select('* from '+table)
    def __call__(self,str,params=None) :
        if len(string.split(str,None,1))==1 :
            list=self.gettab(str)
        else :
            list=self.select(str,params)
        return rdbtab(list, self.description, self.lastsql)
