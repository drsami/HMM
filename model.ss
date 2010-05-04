#lang scheme

(require xml)
(require scheme/list)
(require "data.ss")
(require "xml-functional.ss")

(define lib (make-immutable-hash 
 `((top (,v ,d ,j)))))

(define (makeHMM library)
  (hash-ref library 'top))

(define (makeMonomorphic id emission target)
  `(State ((id ,(number->string id)))
    (transitions ((monomorphic ,(number->string target))))
    (emissions ((monomorphic ,emission)))))

(define (makeIdentity id)
  `(State ((id ,(number->string id)))
          (transitions ((monomorphic ,(number->string id))))
          ,(makeEmissions nucleotides)))

(define (makeIndexedStates id end seqs)
  (foldl (lambda (n acc)
           (let ((res (makeIndexedState (car acc) (car n) (cadr n) makeSubstructure end)))
             (cons (car res) (append (cdr res) (cdr acc)))))
   (list id) seqs))

(define (makeIndexedState id label str substructure terminal)
  (let ((sub (substructure (+ id 1))))
    (let ((nextId (car sub))
          (substates (cdr sub)))
      (cons nextId
      (cons
       `(State ((id ,(number->string id))
                (type "indexed")
                (label ,label)
                (prob "0.999"))
               (emissions ((str ,str)))
               (internalTransition ()
                                   (emission ((nval ,(number->string (+ id 1))) 
                                              (prob "0.0001")))
                                   (emission ((nval ,(number->string (+ id 4))) 
                                              (prob "0.0001")))
                                   (emission ((nval ,(number->string id))
                                              (prob "0.9998"))))
               (terminalTransition ((prob "1.0") (monomorphic ,(number->string terminal)))))
       substates)))))

(define nucleotides '(("A" "0.25") ("T" "0.25") ("C" "0.25") ("G" "0.25")))

(define (makeSubstructure id)
  (list (+ id 6)
   `(State ((id ,(number->string id))
            (NoReset "true"))
           ,(makeEmissions nucleotides)
           ,(makeTransitions `((,(number->string (+ id 1)) "0.99")
                               (,(number->string (- id 1)) "0.01"))))
                            
   `(State ((id ,(number->string (+ id 1)))
            (NoReset "true"))
           ,(makeEmissions nucleotides)
           ,(makeTransitions `((,(number->string (+ id 2)) "0.99")
                               (,(number->string (- id 1)) "0.01"))))
   `(State ((id ,(number->string (+ id 2)))
            (NoReset "true"))
           ,(makeEmissions nucleotides)
           ,(makeTransitions `((,(number->string (- id 1)) "0.99")
                               (,(number->string id) "0.01"))))
   `(State ((id ,(number->string (+ id 3))) 
            (type "silent")
            (NoReset "true")
            (Increment "true"))
           ,(makeTransitions `((,(number->string (+ id 4)) "0.99")
                            (,(number->string id) "0.01"))))
   `(State ((id ,(number->string (+ id 4))) 
            (type "silent")
            (NoReset "true")
            (Increment "true"))
           ,(makeTransitions `((,(number->string (+ id 5)) "0.99")
                            (,(number->string id) "0.01"))))
   `(State ((id ,(number->string (+ id 5)))
            (type "silent")
            (NoReset "true")
            (Increment "true"))
           ,(makeTransitions `((,(number->string id) "0.99")
                            (,(number->string (+ id 3)) "0.01"))))
   ))

(define (makeTransitions em)
  (makeTyped 'transitions 'nval em))

(define (makeEmissions em)
  (makeTyped 'emissions 'val em))

(define (makeTyped type subtype em)
  (append `(,type ())
          (map (lambda (x)
                 `(emission ((,subtype ,(car x))
                             (prob ,(cadr x))))) em)))

(define (makeSilent id label target)
  `(State ((id ,(number->string id))
           ("type" "silent"))
          (transitions ((monomorphic ,(number->string target))))))

(define (list->monomorphic id exit st acc)
  (match st
    [(? empty? st) (list id acc)]
    [_ 
     (let ((target (if (empty? (cdr st)) exit (+ id 1))))
      
     (list->monomorphic
      (+ id 1)
      exit
      (cdr st)
      (append 
        (if (> (string-length (car st)) 1)
            (makeSilent id (car st) target)
            (makeMonomorphic id (car st) target)) acc)))]))

(define (makeGeneral id target dat)
  `(State ((id ,(number->string id))
           (label ,(car dat))
           (type "general"))
          (content () ,(cadr dat))
          (transitions ((monomorphic ,(number->string target))))))
          
(define v-hmm
  (makeIndexedStates 4 1 v))

(define d-hmm
  (makeIndexedStates (car v-hmm) 2 d))

(define j-hmm
  (makeIndexedStates (car d-hmm) 3 j))

(define end-state
  (makeIdentity 3))

(define v-states (map (lambda (x) (xexpr-attr 'id x)) (filter (((curry xexpr-attr?) 'type) "indexed") (cdr v-hmm))))


(define d-states (map (lambda (x) (xexpr-attr 'id x)) (filter (((curry xexpr-attr?) 'type) "indexed") (cdr d-hmm))))

(define j-states (map (lambda (x) (xexpr-attr 'id x)) (filter (((curry xexpr-attr?) 'type) "indexed") (cdr j-hmm))))

(define (makeUnbiasedSilentState id n)
  (let ((numStates (length n)))
    (let ((transitions (map (lambda (x) (list x (number->string (/ 1.0 numStates)))) n)))
      `(State
        ((type "silent")
         (id ,(number->string id)))
        ,(makeTransitions transitions)))))

(define start
  (makeUnbiasedSilentState 0 v-states))

(define v-to-d
  (makeUnbiasedSilentState 1 d-states))

(define d-to-j
  (makeUnbiasedSilentState 2 j-states))

(define hmm
  (concatenate 
   (list `(HMM ()
               ,start
               ,v-to-d
               ,d-to-j
               ,end-state)
         (cdr v-hmm)
         (cdr d-hmm)
         (cdr j-hmm))))    