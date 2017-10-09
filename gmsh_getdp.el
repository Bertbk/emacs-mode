
(require 'cl)

;; Remove previous associations of file extension .geo with idlwave:
(setq auto-mode-alist (delete-if '(lambda (x) (equal (car x) "\\.geo\\'")) auto-mode-alist))

(defvar gmsh/getdp-functions-list
  '("Tanh" "Tan" "Sinh" "Sin" "Sqrt" "Rand" "Modulo" "Log" "Hypot" "Floor" "Fmod" "Fabs" "Exp" "Cosh" "Cos" "Ceil" "Atan" "Atan2" "Asin" "Acos"
    "Printf" "Print" "Sprintf" "Str" "StrCat" "StrFind")
  "Function identifier common to gmsh and getdp")

(defvar gmsh/getdp-keywords-list
  '("For" "EndFor" "If" "EndIf" "Else" "ElseIf" "Include" 
    "Choices" "Label" "Path" "Visible" "Highlight"
    "ReadOnly" "ReadOnlyRange" "Color" "Visible")
  "Keyword identifiers common to gmsh and getdp")

(defvar gmsh/getdp-constants-list
  '("Pi" "View" ".Attributes" ".Light" ".NbIso" ".IntervalsType" ".ShowScale" ".LineWidth" ".ColorTable" "General.AxesFormatX" "General.AxesFormatY" "General.AxesFormatZ" "General.AxesLabelX" "General.AxesLabelY" "General.AxesLabelZ" "General.BackgroundImageFileName" "General.GraphicsFont" "General.GraphicsFontTitle" "General.Axes" "Geometry.CopyMeshingMethod" "Geometry.ExactExtrusion" "Geometry.ExtrudeReturnLateralEntities" "Geometry.ExtrudeSplinePoints" "Geometry.HideCompounds" "Geometry.HighlightOrphans" "Geometry.LabelType" "Geometry.Light" "Geometry.LightTwoSide" "Geometry.Lines" "Geometry.LineNumbers" "Geometry.LineSelectWidth" "Geometry.LineType" "Geometry.LineWidth" "Geometry.MatchGeomAndMesh" "Geometry.Normals" "Geometry.NumSubEdges" "Geometry.OldNewReg" "Geometry.Points" "Geometry.PointNumbers" "Geometry.PointSelectSize" "Geometry.PointSize" "Geometry.PointType" "Geometry.ScalingFactor" "Geometry.OrientedPhysicals" "Geometry.Surfaces" "Geometry.SurfaceNumbers" "Geometry.SurfaceType" "Geometry.Tangents" "Geometry.Tolerance" "Geometry.Transform" "Geometry.TransformXX" "Geometry.TransformXY" "Geometry.TransformXZ" "Geometry.TransformYX" "Geometry.TransformYY" "Geometry.TransformYZ" "Geometry.TransformZX" "Geometry.TransformZY" "Geometry.TransformZZ" "Geometry.Volumes" "Geometry.VolumeNumbers" "Geometry.Color" "Mesh.Algorithm" "Mesh.Algorithm3D" "Mesh.AngleSmoothNormals" "Mesh.AnisoMax" "Mesh.AllowSwapAngle" "Mesh.BdfFieldFormat" "Mesh.Binary" "Mesh.Bunin" "Mesh.Lloyd" "Mesh.SmoothCrossField" "Mesh.CharacteristicLengthFactor""Mesh.CgnsImportOrder" "Mesh.Clip" "Mesh.ColorCarousel" "Mesh.CpuTime" "Mesh.DrawSkinOnly" "Mesh.Dual" "Mesh.ElementOrder" "Mesh.Explode" "Mesh.Format" "Mesh.Hexahedra" "Mesh.HighOrderNumLayers" "Mesh.HighOrderOptimize" "Mesh.HighOrderPoissonRatio" "Mesh.HighOrderSmoothingThreshold" "Mesh.LabelSampling" "Mesh.LabelType" "Mesh.LcIntegrationPrecision" "Mesh.Light" "Mesh.LightLines" "Mesh.LightTwoSide" "Mesh.Lines" "Mesh.LineNumbers" "Mesh.LineWidth" "Mesh.MeshOnlyVisible" "Mesh.MetisAlgorithm" "Mesh.MetisEdgeMatching" "Mesh.MetisRefinementAlgorithm" "Mesh.MinimumCirclePoints" "Mesh.MinimumCurvePoints" "Mesh.MshFileVersion" "Mesh.MshFilePartitioned" "Mesh.MultiplePassesMeshes" "Mesh.PartitionHexWeight" "Mesh.PartitionPrismWeight" "Mesh.PartitionPyramidWeight" "Mesh.PartitionQuadWeight" "Mesh.PartitionTetWeight" "Mesh.PartitionTriWeight" "Mesh.NbHexahedra" "Mesh.NbNodes" "Mesh.NbPartitions" "Mesh.NbPrisms" "Mesh.NbPyramids" "Mesh.NbQuadrangles" "Mesh.NbTetrahedra" "Mesh.NbTriangles" "Mesh.Normals" "Mesh.NumSubEdges" "Mesh.Optimize" "Mesh.OptimizeNetgen" "Mesh.Partitioner" "Mesh.Points" "Mesh.PointNumbers" "Mesh.PointSize" "Mesh.PointType" "Mesh.Prisms" "Mesh.Pyramids" "Mesh.Quadrangles" "Mesh.QualityInf" "Mesh.QualitySup" "Mesh.QualityType" "Mesh.RadiusInf" "Mesh.RadiusSup" "Mesh.RandomFactor" "Mesh.IgnorePartitionBoundary" "Mesh.RecombinationAlgorithm" "Mesh.RecombineAll" "Mesh.Recombine3DAll" "Mesh.RemeshAlgorithm" "Mesh.RemeshParametrization" "Mesh.RefineSteps" "Mesh.Remove4Triangles" "Mesh.ReverseAllNormals" "Mesh.SaveAll" "Mesh.SaveElementTagType" "Mesh.SaveParametric" "Mesh.SaveGroupsOfNodes" "Mesh.ScalingFactor" "Mesh.SecondOrderExperimental" "Mesh.SecondOrderIncomplete" "Mesh.SecondOrderLinear" "Mesh.Smoothing" "Mesh.SmoothNormals" "Mesh.SmoothRatio" "Mesh.SubdivisionAlgorithm" "Mesh.SurfaceEdges" "Mesh.SurfaceFaces" "Mesh.SurfaceNumbers" "Mesh.SwitchElementTags" "Mesh.Tangents" "Mesh.Tetrahedra" "Mesh.ToleranceEdgeLength" "Mesh.Triangles" "Mesh.VolumeEdges" "Mesh.VolumeFaces" "Mesh.VolumeNumbers" "Mesh.Voronoi" "Mesh.ZoneDefinition" "PostProcessing.NbViews" "Mesh.CharacteristicLengthMin" "Mesh.CharacteristicLengthMax")
  "Constant identifiers common to gmsh and getdp")

(defvar gmsh-keywords-list
  '("Return" "Call" "Exit"
    "Include" "SetName" "NonBlockingSystemCall" "SystemCall" "Sleep" "Mesh" "Delete" "BoundingBox" "Merge" "Error"
    "Point" "Physical" "Line" "Loop" "Surface" "Volume" "Compound" "Plane" "Ruled" "Circle" "Ellipse"
    "T2" "T3" "SP" "VP" "TP" "SL" "VL" "TL" "ST" "VT" "TT" "SQ" "VQ" "TQ"
    "SS" "VS" "TS" "SH" "VH" "TH"  "SI" "VI" "TI" "SY" "VY" "TY" "DefineConstant"
    "SetFactory" "Spline" "BSpline" "Block" "Rectangle" "Disk" "Cone" "Wedge" "Sphere" "Box" "Cylinder" "Torus" "ThruSections" "Ruled ThruSections"
    "Compound Surface" "Compound Volume" "Using Wire" "Wire" "Fillet")
  "Gmsh key words")
(defvar gmsh-keywords2-list
  '("newp" "newl" "news" "newv" "newll" "newsl" "newreg"
    "Red" "Blue" "Yellow" "Green" "Orange" "Cyan" "Magneta" "ForestGreen" "SpringGreen" "Gold" "PaleGoldenrod"
    "Pink" "NavyBlue" "Royal Blue" "LightYellow" "LightGreen" "AliceBlue" "Magenta" "Violet" "Orchid" "LightGrey")
  "Gmsh key words2")

(defvar gmsh-functions-list
  '(
    "Extrude" "Dilate" "Rotate" "Symmetry" "Translate" "Boundary" "CombinedBoundary" "Duplicata"
    "Hide" "Show" "TextAttributes"  "Plugin" "Alias" "Combine"
    "BooleanIntersection" "BooleanUnion" "BooleanFragments" "BooleanDifference" "PointsOf"
    )
  "List of function identifiers specific to Gmsh (see also `gmsh/getdp-functions-list').")

(defvar gmsh-constants-list
  '()
  "Gmsh constants.")

(defvar getdp-keywords-list
  '(
    "Type" "TimeFunction" "Value" "Galerkin" "ChangeOfCoordinates" "ChangeOfValues" "InitSolution" "Generate" "GenerateJac" "GenerateOnly" "GenerateOnlyJac" "GenerateSeparate" "SaveSolution" "SaveSolutions" "SetFrequency" "SetTime" "Solve" "SolveJac" "EigenSolve" "IterativeLoop" "IterativeLoopN" "TimeLoopAdaptive" "TimeLoopNewmark" "TimeLoopTheta"
    ;; keywords associated to objects:
    "Analytic" "BasisFunction" "Case" "Constraint" "Criterion" "DefineFunction" "DefineVariable" "DefineConstant" "DefineGroup" "DestinationSystem" "Entity" "EntitySubType" "EntityType" "Equation" "File" "Format" "Formulation" "Frequency" "Function" "FunctionSpace" "GeoElement" "GlobalEquation" "GlobalQuantity" "GlobalTerm" "Group" "In" "IndexOfSystem" "Integral" "Integration" "Jacobian" "Loop" "Name" "NameOfBasisFunction" "NameOfCoef" "NameOfConstraint" "NameOfFormulation" "NameOfMesh" "NameOfPostProcessing" "NameOfSpace" "NameOfSystem" "Node" "NumberOfPoints" "Operation" "OriginSystem" "PostOperation" "PostProcessing" "Quantity" "Region" "RegionRef" "Resolution" "Solver" "SubRegion" "SubRegionRef" "SubSpace" "Support" "Symmetry" "System" "LastTimeStepOnly" "AppendTimeStepToFileName" "OverrideTimeStepValue" "UsingPost" "Append" "MovingBand2D" "MeshMovingBand2D" "InitMovingBand2D" "DeleteFile" "Coefficient" "SendToServer" "Store" "ComputeCommand" "ResolutionChoices" "PostOperationChoices"
    )
  "List of keyword identifiers specific to getdp (see also `gmsh/getdp-keywords-list')")

(defvar getdp-keywords-list-types
  '(
    "All" "Point" "Line" "Triangle" "Quadrangle" "Prism" "Pyramid" "Tetrahedron" "Hexahedron" "Not" 
    "Vol" "VolAxi" "VolAxiRectShell" "VolAxiSphShell" "VolAxiSqu" "VolAxiSquRectShell" "VolAxiSquSphShell" "VolRectShell" "VolSphShell" "Sur" "SurAxi" "Lin" "Form0" "Form1" "Form1P" "Form2" "Form2P" "Form3" 

    ;; Types of:
    "Adapt" "Adaptation" "AliasOf" "Assign" "AssignFromResolution" "AssociatedWith" "BF_CurlEdge" "BF_CurlGroupOfEdges" "BF_CurlGroupOfPerpendicularEdge" "BF_CurlPerpendicularEdge" "BF_dGlobal" "BF_DivFacet" "BF_DivPerpendicularFacet" "BF_Edge" "BF_Facet" "BF_Global" "BF_GradGroupOfNodes" "BF_GradNode" "BF_GroupOfEdges" "BF_GroupOfNodes" "BF_GroupOfPerpendicularEdge" "BF_Node" "BF_NodeX" "BF_NodeY" "BF_NodeZ" "BF_Node_2E" "BF_NodeX_2E" "BF_NodeY_2E" "BF_NodeZ_2E" "BF_One" "BF_PerpendicularEdge" "BF_PerpendicularFacet" "BF_Region" "BF_RegionX" "BF_RegionY" "BF_RegionZ" "BF_Volume" "BF_VolumeX" "BF_VolumeY" "BF_VolumeZ" "BF_Zero" "Break" "DecomposeInSimplex" "Depth" "deRham" "Dimension" "DualEdgesOf" "DualFacetsOf" "DualNodesOf" "DualVolumesOf" "EdgesOf" "EdgesOfTreeIn"  "EigenvalueLegend" "ElementsOf" "Evaluate" "FacetsOf" "FacetsOfTreeIn" "FemEquation" "FourierTransform" "Frequency" "FrequencyLegend"  "Gauss" "GaussLegendre"  "Global" "Gmsh" "GmshParsed" "Gnuplot" "GroupsOfEdgesOf" "GroupsOfEdgesOnNodesOf" "GroupsOfNodesOf" "HarmonicToTime" "If" "Init" "InitFromResolution"  "Integral" "Integral" "Iso" "Lanczos" "Link" "LinkCplx" "Local" "Local" "Network" "NodesOf" "NodeTable" "NoNewLine" "OnBox" "OnElementsOf" "OnGlobal" "OnGrid" "OnLine" "OnPlane" "OnPoint" "OnRegion" "OnSection" "Point" "PostOperation" "Scalar" "SimpleTable" "Skin" "Smoothing"  "Sort" "StartingOn" "StoreInField" "StoreInRegister" "StoreInVariable" "SystemCommand" "Table" "Target"  "TimeLegend" "TimeStep" "TimeTable" "TransferInitSolution" "TransferSolution"  "Update" "UpdateConstraint" "Term" "VolumesOf"  "Local" 
    )
  "List of keyword identifiers specific to getdp (see also `gmsh/getdp-keywords-list')")

(defvar getdp-functions-list
  '("Log10" "Atan2" "Transpose" "TTrace" "Unit" "Complex" "CompX" "CompXX" "CompXY" "CompXZ" "CompY" "CompYX" "CompYY" "CompYZ" "CompZ" "CompZX" "CompZY" "CompZZ" "Im" "Re" "Conj" "Tensor" "TensorDiag" "TensorSym" "TensorV" "Vector" "SquDyadicProduct"
    ;;
    "Dt" "DtDof" "DtDofJacNL" "DtDt" "DtDtDof" "JacNL" "NeverDt" "List" "ListAlt" "h_Jiles" "dhdb_Jiles" "b_Jiles" "dbdh_Jiles"
    ;;
    "Cross" "F_Cos_wt_p" "F_Period" "F_Sin_wt_p" "Fabs" "Fmod" "Hypot" "Norm" "SquNorm" "SurfaceArea"
    ;; fields:
    "BF" "Curl" "CurlInv" "dInv" "Div" "DivInv" "Dof" "Grad" "GradInv" "Rot" "RotInv"
    ;; current Values:
    ;; "$DTime" "$EigenvalueImag" "$EigenvalueReal" "$Iteration" "$Theta" "$Time" "$TimeStep" "$X" "$XS" "$Y" "$YS" "$Z" "$ZS"
    ;; "X[]" "Y[]" "Z[]" "XYZ[]"
    ;; misc functions:
    "dInterpolationAkima" "dInterpolationLinear" "F_CompElementNum" "Order" "Rand"
    ;; Green functions:
    "Helmholtz" "Laplace"
    )
  "List with function identifiers specific to getdp (see also `gmsh/getdp-functions-list').")

(defvar getdp-constants-list
  '("0D" "1D" "2D" "3D" 
    ;; current Values:
    "$DTime" "$EigenvalueImag" "$EigenvalueReal" "$Iteration" "$Theta" "$Time" "$TimeStep" "$X" "$XS" "$Y" "$YS" "$Z" "$ZS" "X[]" "Y[]" "Z[]" "XYZ[]" 
    )
  "Getdp constants and current values.")

(defvar gmsh/getdp-block-statements '(("For" "EndFor")
				      ("If" "EndIf")
                                      ("If" "Else")
                                      ("If" "ElseIf")
                                      ("Else" "EndIf")
                                      ("ElseIf" "ElseIf")
                                      ("ElseIf" "Else")
                                      ("ElseIf" "EndIf")
                                      )
  "List of block statements that are common to gmsh and getdp.
This list is mainly used for indentation.  Each entry is a list
with the begin and end statement for a block.")

(defvar gmsh-block-statements  '(("Function" "Return"))
  "Gmsh block statements. Mainly for indentation.
In each entry the first element is the block start the second the block end.")

(defvar getdp-block-statements '(("Macro" "Return"))
  "GetDP block statements. Mainly for indentation.
In each entry the first element is the block start the second the block end.")

(defun gmsh/getdp-highlight-In-after-For (lim)
  "Find In after For"
  (interactive)
  (let (found)
    (catch 0
      (while (setq found (search-forward-regexp "\\_<In\\_>" lim 'noErr))
	(and
	 (null (syntax-ppss-context (syntax-ppss)))
	 (let ((inhibit-changing-match-data t)) (looking-back "\\_<For[[:blank:]]+\\(\\sw\\|\\s_\\)+[[:blank:]]+In"))
	 (throw 0 found))))))

(defun getdp-highlight-Vector (lim)
  "Find Vector after Type"
  (interactive)
  (let (found)
    (catch 0
      (while (setq found (search-forward-regexp "\\_<Vector\\_>" lim 'noErr))
	(and
	 (null (syntax-ppss-context (syntax-ppss)))
	 (let ((inhibit-changing-match-data t)) (looking-back "\\_<Type[[:blank:]]+Vector"))
	 (throw 0 found))))))

(defun getdp-highlight-Region (lim)
  "Find Region after EntityType"
  (interactive)
  (let (found)
    (catch 0
      (while (setq found (search-forward-regexp "\\_<Region\\_>" lim 'noErr))
	(and
	 (null (syntax-ppss-context (syntax-ppss)))
	 (let ((inhibit-changing-match-data t)) (looking-back "\\_<EntityType[[:blank:]]+Region"))
	 (throw 0 found))))))




(defvar indent-amount 2
  "Number of columns to insert additionally for each indentation level.")

(defun nonspace-before (p)
  "Return list (c p l) with character c position p and number of
lines l when traveling backwards stopping at position p behind
the first non-space character c. Thereby, l is the number of
lines passed on the way."
  (let (c (l 0))
    (while (and
	    (setq c (char-before p))
	    (or (when (= c ?\n) (setq l (1+ l)))
		(= (char-syntax c) ?\ )))
      (setq p (1- p)))
    (list c p l)))

(defvar indent-relative-blocks nil
  "List of block statement descriptions recognized by `indent-relative-function'.
Each block statement description is itself a list containing two
strings.  The first string is the identifier for the beginning of
the block the second that one for the end.")

(make-variable-buffer-local 'indent-relative-blocks)

(defun indent-relative-function ()
  "Indent current line relative to previous one."
  (interactive)
  (save-excursion
    (if (eq (syntax-ppss-context (syntax-ppss ;; note: syntax-ppss can move point!
				  (line-beginning-position))) 'string)
	'noindent ;; line starts in string ==> return 'noindent
      ;; calculate indent level
      ;; search for the previous non-empty line:
      (let (b e prevIndent relDepth (p (point))
	      to
	      (re-block (regexp-opt (apply 'append indent-relative-blocks) 'words))
	      (re-end (concat "\\(\\s)\\|" (regexp-opt (mapcar 'cadr indent-relative-blocks) 'words) "\\)")))
	(save-excursion
	  (setq relDepth (syntax-ppss-depth (syntax-ppss)))
	  (beginning-of-line (- 1 (nth 2 (nonspace-before (line-beginning-position)))))
	  (setq prevIndent (current-indentation))
	  (move-to-column prevIndent)
	  (setq relDepth (- relDepth (syntax-ppss-depth (syntax-ppss))))
	  (when (and
		 (null (member (syntax-ppss-context (syntax-ppss)) '(string comment)))
		 (looking-at-p re-end))
	    (setq relDepth (1+ relDepth)))
	  ;; look for block keywords:
	  (setq e (line-end-position))
	  (while (search-forward-regexp re-block e 'noErr)
	    (unless (syntax-ppss-context (syntax-ppss))
	      (if (assoc (match-string-no-properties 0) indent-relative-blocks)
		  (setq relDepth (1+ relDepth))
		(setq relDepth (1- relDepth))
		))))
	(goto-char p)
	(beginning-of-line)
	(delete-horizontal-space)
	(when (looking-at-p re-end)
	  (setq relDepth (1- relDepth)))
	(setq to (+ prevIndent (* relDepth indent-amount)))
	(unless (= to (current-indentation))
	  (indent-to to)))))
  (when (< (current-column) (current-indentation))
    (move-to-column (current-indentation))))

(defmacro key-with-indent (k)
  "Define k as key with indent."
  (let ((sym (make-symbol (concat "electric-char-" (char-to-string k)))))
  `(progn
     (defun ,sym ()
       "Run `self-insert-command' and `indent-relative-function' for this character."
       (interactive)
       (self-insert-command 1)
       (indent-relative-function)
       )
     (local-set-key [,k] (quote ,sym))
     )))

(defun gmsh/getdp-common-settings ()
  "Common settings for gmsh and getdp."
  (setq font-lock-defaults (list
			    ;;;;;;;;;;
			    ;; KEYWORDS:
			    (list
			     ;;(cons (regexp-opt gmsh/getdp-keywords-list 'symbols) 'font-lock-keyword-face) ;; purple
                             (cons (regexp-opt gmsh/getdp-keywords-list 'symbols) 'font-lock-function-name-face) ;; purple
			     (cons (regexp-opt gmsh/getdp-functions-list 'symbols) 'font-lock-function-name-face) ;; blue
			     (cons (regexp-opt gmsh/getdp-constants-list 'symbols) 'font-lock-constant-face)
			     ;;(cons 'gmsh/getdp-highlight-In-after-For 'font-lock-keyword-face)
                             (cons 'gmsh/getdp-highlight-In-after-For 'font-lock-function-name-face)
                             (cons 'getdp-highlight-Vector 'font-lock-type-face)
                             (cons 'getdp-highlight-Region 'font-lock-type-face))
			    ;;;;;;;;;;
			    nil ; keywords-only
			    nil ; case-fold
			    ;;;;;;;;;;
			    ;; syntax-list:
			    '(
			      ("_." . "w")
			      (?/ . ". 1456")
			      (?/ . ". 124b")
			      (?* . ". 23")
			      (?\n . "> b")
			      (?= . ".") (?+ . ".") (?- . ".")
			      (?$ . "_")
			      )
			    ))
  (local-set-key (kbd "<return>") '(lambda () (interactive) (newline) (indent-relative-function)))
  (key-with-indent ?}) (key-with-indent ?{)
  (key-with-indent ?() (key-with-indent ?))
  (key-with-indent ?[)  (key-with-indent ?])
  (setq indent-line-function 'indent-relative-function))

(if (version-list-< (version-to-list emacs-version) '(24 0))
    (defadvice regexp-opt (after symbols activate)
      (if (eq paren 'symbols)
	  (setq ad-return-value (concat "\\_<" ad-return-value "\\_>")))))

(defun last-cons (list)
  "Return last cons in list. Note, that list must contain at least one element."
  (unless (and list (listp list))
    (error "Call of last-cons with non-list or empty list argument."))
  (while (cdr list)
    (setq list (cdr list)))
  list)

(define-derived-mode gmsh-mode c++-mode "gmsh"
  "Major mode for editing gmsh geometry definitions."
  (gmsh/getdp-common-settings)
  (setcdr (last-cons (car font-lock-defaults))
	  (list (cons (regexp-opt gmsh-keywords-list 'symbols) 'font-lock-keyword-face) ;; purple
                (cons (regexp-opt gmsh-keywords2-list 'symbols) 'font-lock-type-face) ;; green
		(cons (regexp-opt gmsh-functions-list 'symbols) 'font-lock-function-name-face) ;; blue
		(cons (regexp-opt gmsh-constants-list 'symbols) 'font-lock-constant-face)
		))
  (setq indent-relative-blocks (append gmsh/getdp-block-statements gmsh-block-statements))
  (font-lock-mode 1)
  )

(define-derived-mode getdp-mode c++-mode "getdp"
  "Major mode for editing getdp geometry definitions."
  (gmsh/getdp-common-settings)
  (setcdr (last-cons (car font-lock-defaults))
	  (list (cons (regexp-opt getdp-keywords-list 'symbols) 'font-lock-keyword-face)
                (cons (regexp-opt getdp-keywords-list-types 'symbols) 'font-lock-type-face)
                (cons (regexp-opt getdp-functions-list 'symbols) 'font-lock-function-name-face) ;; blue
		(cons (regexp-opt getdp-constants-list 'symbols) 'font-lock-constant-face)
		))
  (setq indent-relative-blocks (append gmsh/getdp-block-statements getdp-block-statements))
  )

(add-to-list 'auto-mode-alist '("\\.pro\\'" . getdp-mode))
(add-to-list 'auto-mode-alist '("\\.geo\\'" . gmsh-mode))

(require 'info-look)

(info-lookup-add-help
 :mode 'getdp-mode
 :regexp "[[:alnum:]_]+"
 :doc-spec '(("(getdp)Syntax Index" nil "")))

(info-lookup-add-help
 :mode 'gmsh-mode
 :regexp "[[:alnum:]_]+"
 :doc-spec '(("(gmsh)Syntax Index" nil "")))

(defmacro stl2geo-vertex (pts cnt str name scale)
  "Read one vertex from stl file. Put it into hashtable pts with
index value cnt if it does not yet exist there else return its
index value in pts. Point moves behind the vertex."
  `(progn
     (search-forward "vertex")
     (let* ((b (current-buffer))
	    (pt (list (read b) (read b) (read b)))
	    (idx (gethash pt ,pts)))
       (unless idx
	 (puthash pt ,cnt ,pts)
	 (setq idx ,cnt)
	 (setq ,str (concat ,str "\np" ,name "[" (number-to-string idx) "]=newp; Point(p" ,name "[" (number-to-string idx) "])={"
			    (number-to-string (* ,scale (nth 0 pt))) ","
			    (number-to-string (* ,scale (nth 1 pt))) ","
			    (number-to-string (* ,scale (nth 2 pt))) ",mshSize};"))
	 (setq ,cnt (1+ ,cnt)))
       idx)))

(defmacro stl2geo-edge (edges p1 p2 cnt str name)
  `(let* (idx)
     (setq idx (gethash (cons ,p1 ,p2) ,edges))
     (unless idx
       (setq idx (gethash (cons ,p2 ,p1) ,edges))
       (when idx (setq idx (- idx))))
     (unless idx
       (puthash (cons ,p1 ,p2) ,cnt ,edges)
       (setq idx ,cnt)
       (setq ,str (concat ,str "\nl" ,name "[" (number-to-string idx) "]=newl; Line(l" ,name "[" (number-to-string idx)
			  "])={p" ,name "[" (number-to-string ,p1)
			  "],p" ,name "[" (number-to-string ,p2) "]};"))
       (setq ,cnt (1+ ,cnt)))
     idx))

(defun stl2geo (name scale geo)
  "Transformes stl grid of current buffer into gmsh geo file."
  (interactive "sName of entity:\nnScale factor:\nBGeo-buffer (without extension .geo):")
  (save-excursion
    (goto-char (point-min))
    (let* (facet
	   str
	   lineloop
	   (cntPts 0)
	   (cntEdges 1)
	   (cntLineLoop 0)
	   (size (/ (line-number-at-pos (point-max)) 7))
	   (pts (make-hash-table :test 'equal :size size))
	   (edges (make-hash-table :test 'equal :size size))
	   (geobuf (get-buffer-create geo)))
      (while (search-forward-regexp "^outer loop" nil 'noErr)
	(setq str nil)
	(setq facet (list (stl2geo-vertex pts cntPts str name scale)
			  (stl2geo-vertex pts cntPts str name scale)
			  (stl2geo-vertex pts cntPts str name scale)))
	(setq lineloop (list (stl2geo-edge edges (nth 0 facet) (nth 1 facet) cntEdges str name)
			     (stl2geo-edge edges (nth 1 facet) (nth 2 facet) cntEdges str name)
			     (stl2geo-edge edges (nth 2 facet) (nth 0 facet) cntEdges str name)))
	(setq str (concat str
			  "\nll" name "[" (number-to-string cntLineLoop) "]=newll; Line Loop(ll" name "[" (number-to-string cntLineLoop) "])={"
			  (if (< (nth 0 lineloop) 0 ) "-l" "l") name "[" (number-to-string (abs (nth 0 lineloop))) "],"
			  (if (< (nth 1 lineloop) 0 ) "-l" "l") name "[" (number-to-string (abs (nth 1 lineloop))) "],"
			  (if (< (nth 2 lineloop) 0 ) "-l" "l")  name "[" (number-to-string (abs (nth 2 lineloop))) "]};\n"
			  "s" name "[" (number-to-string cntLineLoop) "]=news; Plane Surface(s" name "[" (number-to-string cntLineLoop) "])={ll" name "[" (number-to-string cntLineLoop) "]};"
			  ))
	(setq cntLineLoop (1+ cntLineLoop))
	(with-current-buffer geobuf
	  (insert str)
	  )))))

(defun gmsh-import-stl (stl-file-name geo-entity-name scale)
  "Insert ASCII-STL mesh from file `stl-file-name' at point.
The geometrical entities are suffixed by `geo-entity-name'. The points are scaled by `scale'."
  (interactive "fSTL file:
sName suffix for geometrical entities:
nScale factor:")
  (let ((geo-buf (current-buffer)))
    (with-temp-buffer
      (insert-file-contents stl-file-name)
      (stl2geo geo-entity-name scale geo-buf))))

(easy-menu-define nil gmsh-mode-map
  "Menu for gmsh (import of ascii-stl)."
  '("Gmsh"
    ["Import ASCII STL" gmsh-import-stl t]))

(provide 'gmsh)

