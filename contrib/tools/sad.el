(defvar sad-mode-map (make-sparse-keymap))
(define-key sad-mode-map "\t" 'sad-indent-line)

(define-generic-mode 'sad-mode
  '("!" ("(*" . "*)"))
  '()
  '(("\\<\\(Class\\|Module\\|Block\\|With\\|Do\\|While\\|For\\|If\\|Switch\\|Which\\|Return\\|Check\\|Catch\\|Throw\\|Break\\|Continue\\|Goto\\|Label\\|Exit\\|On\\|Off\\|Clear\\|SetAttributes\\|Constructor\\|And\\|Or\\|Not\\|Begin\\|End\\|BeginPackage\\|EndPackage\\) *\\["
     (1 font-lock-keyword-face))
    ("\\<\\(Length\\|Thread\\|Table\\|Map\\|Apply\\|Thread\\|Scan\\|MapThread\\|ScanThread\\|Cases\\|Select\\|Position\\|MapAt\\|ReplaceAt\\|Sum\\|Product\\|SwitchCases\\|SelectCases\\|Count\\|Rest\\|First\\|Second\\|Third\\|Take\\|Drop\\|Join\\|Append\\|Prepend\\|Complement\\|Range\\|Union\\|Sort\\|Override\\|AppendTo\\|PrependTo\\|ReplacePart\\|Head\\|Flatten\\|Partition\\|DeleteCases\\|OpenRead\\|OpenWrite\\|OpenAppend\\|Read\\|Write\\|ReadString\\|WriteString\\|Close\\|Get\\|RotateLeft\\|RotateRight\\) *\\["
     (1 font-lock-function-name-face))
    ("\\<\\(Min\\|Max\\|MinMax\\|Restrict\\|Abs\\|Floor\\|Ceiling\\|Round\\|Sqrt\\|Mod\\|Plus\\|Times\\|Dot\\|Det\\|Transpose\\|LinearSolve\\|Eigensystem\\|SingularValues\\|Fourier\\|Complex\\|Sin\\|Cos\\|Tan\\|Sinh\\|Cosh\\|Tanh\\|ArcSin\\|ArcCos\\|ArcTan\\|ArcSinh\\|ArcCosh\\|ArcTanh\\|Log\\|Exp\\|Random\\|GaussRandom\\|Bessel[JYIK]\\|Erf\\|Factorial\\|D\\|Tr\\) *\\["
     (1 font-lock-type-face))
    ("\\<\\(True\\|False\\|Null\\|Undefined\\|EndOfFile\\|Infinity\\|INF\\|NaN\\|Pi\\|Degree\\|E\\|I\\|Automatic\\|None\\|Both\\|ElectronRadius\\|ElectronCharge\\|SpeedOfLight\\|PlanckConstant\\|FineStructureConstant\\|ElectronMass\\|ProtonRadius\\|GoldenRatio\\|Real\\|String\\|Constant\\|Hold\\|Literal\\|HoldPattern\\|HoldFirst\\|HoldAll\\|HoldRest\\|HoldNone\\|Word\\|Expression\\|WordSeparator\\|ReadNewRecord\\|NullWords\\|$Failed\\|FormatType\\|InputForm\\|HoldForm\\|GenericSymbolForm\\|Default\\|EulerGamma\\|SIMu0\\|MKSAMu0\\|SIEpsilon0\\|MKSAEpsilon0\\|BuiltinFunction\\|ProtonMass\\)\\($\\|[^a-zA-Z0-9$]\\)"
     (1 font-lock-constant-face))
    ("\\({\\|}\\|,\\|\\;\\|_\\)"
     (1 font-lock-comment-delimiter-face))
    ("\\(:=\\)"
     (1 font-lock-preprocessor-face))
    ("[^=<>]\\(=\\)[^=><]"
     (1 font-lock-preprocessor-face))
    ("\\(->\\|:>\\)"
     (1 font-lock-preprocessor-face))
    )
  '("\\.sad$" "\\.n$")
  (list `sad-setup)
  "Mode for SAD files.")

(defun sad-setup ()
  (use-local-map sad-mode-map)
  (make-local-variable 'indent-line-function)
  (setq indent-line-function 'sad-indent-line)
  (make-local-variable 'parse-sexp-ignore-comments)
  (setq parse-sexp-ignore-comments t)
  )


(defun sad-indent-line ()
  (interactive)
  (let ((p (- (point-max) (point))))
    (beginning-of-line)
    (forward-to-indentation 0)
    (if (not (string= (buffer-substring (point) (1+ (point))) "!"))
        (let ((p1 (point)))
          (sad-eol-to-nocomment)
          (if (< (point) p1)
              (if (bolp)
                  (sad-indent-up-list p1)
                (let ((c (buffer-substring (1- (point)) (point))))
                  (or (string= c "\\")
                      (if (condition-case nil
                              (or (string-match "[0-9,A-Z,a-z]" c)
                                  (string-match c ")}];,"))
                            (error nil))
                          (sad-indent-up-list p1)
                        (let ((n (+ 2 (back-to-indentation1))))
                          (goto-char p1)
                          (indent-line-to n)))))))))
    (goto-char (max (progn (beginning-of-line) (point))
                    (- (point-max) p)))))

(defun sad-eol-to-nocomment ()
  (let ((q t))
    (while q
      (end-of-line 0)
      (if (= (point) 0) (setq q nil)
        (let ((p (point)))
          (back-to-indentation1)
          (or (eolp)
              (if (not (string= (buffer-substring (point) (1+ (point))) "!"))
                  (setq q nil))))))
    (beginning-of-line)
    (if (not (eolp))
        (re-search-forward "\\( *!.*\\| *\\)$" nil t))))
          
(defun sad-indent-up-list (p1)
  (goto-char p1)
  (condition-case nil
      (progn     
        (backward-up-list 1)
        (if (= p1 (point))
            (sad-indent-out-of-list p1)
          (let ((n (+ 2 (back-to-indentation1))))
            (goto-char p1)
            (indent-line-to n))))
    (error (progn 
             (sad-indent-out-of-list p1)
             nil))))

(defun sad-indent-out-of-list (p1)
  (let ((q t))
    (while q
      (sad-eol-to-nocomment)
      (beginning-of-line)
      (let ((p (point)))
        (condition-case nil (progn (backward-up-list 1)
                               (setq q (< (point) p)))
            (error (progn (setq q nil) nil))))))
  (let ((n (back-to-indentation1)))
    (goto-char p1)
    (indent-line-to n)))

(defun back-to-indentation1 nil
  (let ((p2 (progn (beginning-of-line) (point))))
    (back-to-indentation)
    (- (point) p2)))
