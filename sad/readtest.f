FFS;

f=OpenRead["~/SAD/oldsad/sad/examples.sad"];
Read[f,String]
Read[f,Word]
Read[f,Word*10,ReadNewRecord->False]
Read[f,Word]
Read[f]
ReadString[f]
Close[f];

end





