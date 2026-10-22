;;; csvtab2org --- A generic CSV to Org-table converter -*- lexical-binding: t; -*-

;;; Commentary:

;;; Code:

;;; my-org-table-from-csv
(defun my-org-table-from-csv (csv-file caption
                              &optional header-vector
                                        numeric-cols
                                        first-col-values
                                        floatfmt
                                        text-size)
  "Return an Org‑mode table (as a string) built from CSV-FILE.

CSV-FILE : path to the CSV file.
CAPTION : caption string.
HEADER-VECTOR : optional list of column titles.
NUMERIC-COLS : optional list of zero-based column indices.
FIRST-COL-VALUES : optional list to replace the first column values.
FLOATFMT : numeric format string e.g. \".2f\" (default \".1f\").
TEXT-SIZE : LaTeX size command e.g. \"\\small\" (optional)."

  ;; Check file existence
  (unless (file-exists-p csv-file)
    (error "CSV file %s does not exist" csv-file))

  (let* (;; Numeric format
         (float-format (if (and floatfmt (string-match "^\\.\\([0-9]+\\)f$" floatfmt))
                           (format "%%.%sf" (match-string 1 floatfmt))
                         "%.1f"))

         ;; Read and parse the CSV
         (raw-data
          (with-temp-buffer
            (insert-file-contents csv-file)
            ;; Parse line by line (simple and reliable method)
            (let ((lines '()))
              (goto-char (point-min))
              (while (not (eobp))
                (let ((line (buffer-substring-no-properties 
                             (line-beginning-position) 
                             (line-end-position))))
                  (when (> (length (string-trim line)) 0)
                    ;; Split on commas, handle simple quoted fields
                    (push (split-string line "," t "[ \t\n]+") lines))
                  (forward-line 1)))
              (nreverse lines))))

         ;; Header: first line or header-vector
         (final-header (or header-vector (nth 0 raw-data)))

         ;; Data rows (without the first line)
         (data-rows (if header-vector
                        raw-data
                      (cdr raw-data)))
         (n-rows (length data-rows))

         ;; Verify first-col-values length
         (first-col-values (when first-col-values
                             (let ((lst (if (vectorp first-col-values)
                                            (append first-col-values nil)
                                          first-col-values)))
                               (unless (= (length lst) n-rows)
                                 (error "Length of FIRST-COL-VALUES (%d) does not match number of rows (%d)"
                                        (length lst) n-rows))
                               lst)))

         ;; Numeric columns
         (num-cols (or numeric-cols '()))

         ;; Build formatted rows
         (formatted-rows
          (cl-loop for row in data-rows
                   for i from 0
                   collect
                   (cl-loop for j from 0
                            for cell in row
                            collect
                            (cond
                             ;; Replace first column if requested
                             ((and (= j 0) first-col-values (nth i first-col-values))
                              (nth i first-col-values))
                             ;; Format numeric columns
                             ((member j num-cols)
                              (format float-format (string-to-number cell)))
                             ;; Otherwise keep the value as-is
                             (t cell)))))

         ;; Combine header and data rows
         (all-rows (cons final-header formatted-rows))

         ;; Calculate column widths for alignment
         (n-columns (length final-header))
         (column-widths
          (cl-loop for col from 0 below n-columns
                   collect
                   (cl-loop for row in all-rows
                            for cell = (nth col row)
                            when cell
                            maximize (length (format "%s" cell)))))

         ;; Pad cells to align columns
         (pad-cell (lambda (cell width)
                     (let ((str (format "%s" cell)))
                       (concat str (make-string (- width (length str)) ?\s)))))

         ;; Build aligned rows (using cl-loop for indexing)
         (aligned-rows
          (cl-loop for row in all-rows
                   collect
                   (concat "| "
                           (cl-loop for cell in row
                                    for col from 0
                                    collect (funcall pad-cell cell (nth col column-widths))
                                    into padded-cells
                                    finally (return (mapconcat 'identity padded-cells " | ")))
                           " |")))

         ;; Build separator line
         (separator
          (concat "|-"
                  (mapconcat (lambda (w) (make-string w ?-)) column-widths "-|")
                  "-|"))

         ;; Assemble the final output
         (attr-line (if text-size
                        (format "#+ATTR_LATEX: :booktabs t :font %s\n" text-size)
                      "#+ATTR_LATEX: :booktabs t\n")))

    ;; Return the complete table
    (concat
     (format "#+CAPTION: %s\n" caption)
     attr-line
     (nth 0 aligned-rows) "\n"
     separator "\n"
     (mapconcat 'identity (cdr aligned-rows) "\n")
     "\n")))

;; ## **Exemple d'utilisation :**

;; ```elisp
;; ;; Créer un fichier test
;; (with-temp-file "/tmp/test.csv"
;;   (insert "Nom,Age,Ville\n"
;;           "Alice,30,Paris\n"
;;           "Bob,25,Lyon\n"))

;; ;; Générer le tableau
;; (my-org-table-from-csv 
;;  "/tmp/test.csv"
;;  "Exemple de tableau"
;;  nil  ; header-vector (utilise la 1ère ligne)
;;  '(1) ; colonnes numériques (Age)
;;  nil  ; first-col-values
;;  ".0f" ; format : entier
;;  "\\small") ; taille LaTeX
;; ```

;;; csvtab2org.el ends here
