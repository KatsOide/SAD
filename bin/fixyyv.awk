/(yypvt|yyvsp)\[[^\]]*\]/{
   printf "      "
   for (i=1;i<=NF;i++){
     if($i ~ /(yypvt|yyvsp)\[[^\]]*\]/) {
          p=index($i,"[")+1
          $i=substr($i,1,p-7)"yyv(yypvt"substr($i,p,index($i,"]")-p)")"substr($i,index($i,"]")+1)
     }
     printf " %s", $i
   }
   printf "\n"
   next
}
  { print $0}
