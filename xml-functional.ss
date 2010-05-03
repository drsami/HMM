(module xml-functional scheme
  (provide
   file->xexpr
   xml-strip-headers
   xexpr-filter
   xexpr-map
   xexpr-map/tag
   xexpr-map/attr
   xexpr-map/cdata
   xexpr-is?
   xexpr-is?/c
   xexpr-attr
   xexpr-attr!
   xexpr-attr?
   concatenate)
  (require xml)
  (require scheme/function)
  (require scheme/file)
  (require scheme/system)
  ;(require srfi/13)
  ;(require srfi/1)

;; file->xexpr
;; input:  a path containing xml
;; output: an xexpr representation of the xml
(define (file->xexpr fn)
  (string->xexpr (xml-strip-headers (file->string fn))))
  
;; xml-strip-headers
;; input:  a string, possibly with xml encoding headers
;; output: the string without headers
(define (xml-strip-headers str)
  (let ((pstr (regexp-split #rx"<\\?xml.*\\?>" str)))
    (if (< 1 (length pstr))
        (cadr pstr)
        (car pstr))))

;; xexpr-map
;; input:  f (xexpr -> xexpr), an xexpr
;; output: an xexpr
(define (xexpr-map f xexpr)
  (match xexpr
    [(? string?) (f xexpr)]
    [(list t attr children ...) (append (f xexpr) (map ((curry xexpr-map) f) children))]))

;; xexpr-map/tag: map over the tags of an xexpr
;; input:  (tag -> tag) xexpr
;; output: xexpr
(define (xexpr-map/tag f xexpr)
  (xexpr-map
   (lambda (x)
     (if (list? x)
         (list (f (car x)) (cadr x))
         x))
   xexpr))

;; xexpr-map/attr: map over the attributes of an xexpr
;; input:  (attr -> attr) xexpr
;; output: xexpr
(define (xexpr-map/attr f xexpr)
  (xexpr-map
   (lambda (x)
     (if (list? x)
         (list (car x) (f (cadr x)))
         x))
     xexpr))

;; xexpr-map/cdata: map over the CDATA of an xexpr
;; intput: (string -> string) xexpr
;; output: xexpr
(define (xexpr-map/cdata f xexpr)
  (xexpr-map 
   (lambda (x) 
     (if (string? x)
         (f x)
         (list (car x) (cadr x)))) xexpr))

;; xexpr-is?
;; input:  tag xexpr
;; output: #t if the xexpr is of type tag, otherwise #f
(define (xexpr-is? tag xexpr)
  (if (and (list? xexpr) (not (empty? xexpr)))
      (equal? (car xexpr) tag)
      #f))

;; xexpr-is?/c
;;  a curried version of xexpr-is?
(define (xexpr-is?/c tag)
  (lambda (xexpr)
    (xexpr-is? tag xexpr)))


;; xexpr-attr
;; input:  (attr symbol) xexpr
;; output: the value of the attribute 
(define (xexpr-attr attr xexpr)
  (match xexpr
    [(? string?) (error "CDATA cannot have attributes")]
    [(list t attrList children ...) (extract-attr attr attrList)]
    [(list t attrList) (extract-attr attr attrList)]))

(define (extract-attr attr lst)
  (match lst
    [(list (list k v) rst ...) 
     (if (equal? k attr)
         v
         (extract-attr attr rst))]
    [(list (list k v))
     (if (equal? k attr)
         v
         "")]
    [(list) ""]))
  
;; xexpr-attr!
;; input:
;; output: 
;; I'm not checking to see if an attribute is already set
;; shame on me.
(define (xexpr-attr! attr val xexpr)
  (match xexpr
    [(list tag attrl rst)
     (list tag (cons (list attr val) attrl))]
    [(list tag attrl rst ...)
     (list tag (cons (list attr val) attrl) rst)]))
  
  
(define (xexpr-attr? attr val xexpr)
  (match xexpr
    [(? string?) (error "CDATA cannot have attributes")]
    [(list t attrList children ...) (has-attr? attr val attrList)]
    [(list t attrList) (has-attr? attr val attrList)]))
  
(define (has-attr? attr val lst)
  (match lst
    [(list (list k v) rst ...)
     (if (equal? k attr)
         (equal? v val)
         (has-attr? attr val rst))]
    [(list (list k v))
     (if (equal? k attr)
         (equal? v val)
         #f)]
    [(list) #f]))
  
(define (concatenate lists)
  (foldr append '() lists))
  
(define (xexpr-filter f x)
  (match x
    [(? string?) '()]
    [(list) '()]
    [(list t attr children ...) 
     (if (f x)
         (cons x 
               (concatenate (map (lambda (x) (xexpr-filter f x)) children)))
         (concatenate (map (lambda (x) (xexpr-filter f x)) children)))]))
  )