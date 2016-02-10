c     Don't confuse, Emacs. This is -*- fortran -*- mode!
      logical errflg
      common /codcor/eptsol,errflg(4),errval(4),dpshft,optiv(18),
     $               ibckup,nmona,nmonact,nstra,nstract,
     $               itmon,itemon,itstr,itestr
c          mcepst msolvg mcmon mcnrmc mdpmax palgn pcbak pundo prkick           
c eptsol     o      o                          o                                
c errflg                   o                                                    
c errval                   o     o             o                                
c dpshft                         o      o                                       
c ibckup                                             o      o     o             
c optiv                                              o                          
c nmona                    o     o             o                                
c nmonact                                                                       
c nstra                                        o                  o             
c nstract                                                                       
c itmon                                                                         
c itemon                                                                        
c itstr                                                                         
c itestr                                                                        
c          pbump  corinit mcstr mccor mcrcod mcrmda monel pwrite mrecal         
c eptsol     o      o             o                                             
c errflg            o                                                           
c errval            o                                                           
c dpshft            o                                                           
c ibckup     o      o                                                           
c optiv                                                                         
c nmona      o                    o     o      o      o     o      o            
c nmonact                                                                       
c nstra                     o     o                         o      o            
c nstract                                                                       
c itmon             o                                                           
c itemon            o                                                           
c itstr             o       o                                                   
c itestr            o       o                                                   
c          pkill pstati mstore ptrim petcod mclear twsdrw pbumps monact         
c eptsol                                                                        
c errflg                                                                        
c errval                                                                        
c dpshft                                                                        
c ibckup                                                                        
c optiv                                                                         
c nmona      o     o      o      o      o            o                          
c nmonact                                                          o            
c nstra      o     o      o                   o              o                  
c nstract                                                                       
c itmon                                                            o            
c itemon                                                                        
c itstr                                                                         
c itestr                                                                        
c          mcrmat mbmp msolb pcrmat ptol pcset mstack mweght preadmon
c eptsol                                                                        
c errflg                                                                        
c errval                                                                        
c dpshft                                                                        
c ibckup                                                                        
c optiv                                                         o                
c nmona                                          o      o                       
c nmonact                                                       o
c nstra      o     o     o     o     o     o     o                              
c nstract                                                                       
c itmon                                                                         
c itemon                                                                        
c itstr                                                                         
c itestr                                                                        
c          preadstr pwrtmon pwrtstr
c eptsol                                                                        
c errflg                                                                        
c errval                                                                        
c dpshft                                                                        
c ibckup                                                                        
c optiv       
c nmona               o
c nmonact                                                                       
c nstra      o                o
c nstract                                                                       
c itmon                                                                         
c itemon                                                                        
c itstr                                                                         
c itestr                                                                        
