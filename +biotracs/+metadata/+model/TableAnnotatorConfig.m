%"""
%biotracs.metadata.model.TableAnnotatorConfig
%Configuration object.
%* License: BIOASTER License
%* Created by: Bioinformatics team, Omics Hub, BIOASTER Technology Research Institute (http://www.bioaster.org), 2017
%* See: biotracs.metadata.model.TableAnnotator
%"""

classdef TableAnnotatorConfig < biotracs.core.mvc.model.ProcessConfig
	 
	 properties(Constant)
	 end
	 
	 properties(SetAccess = protected)
	 end

	 % -------------------------------------------------------
	 % Public methods
	 % -------------------------------------------------------
	 
	 methods
		  
		  % Constructor
		  function this = TableAnnotatorConfig( )
				this@biotracs.core.mvc.model.ProcessConfig();
                this.createParam('TokenDelimiter', '_', 'Constraint', biotracs.core.constraint.IsText());
                this.createParam('KeyValueDelimiter', ':', 'Constraint', biotracs.core.constraint.IsText());
                this.createParam('HeadersToUse', {}, 'Constraint', biotracs.core.constraint.IsText('IsScalar', false));
                this.createParam('UserIdHeader', 'UserId');
                this.createParam('AbbreviationTable', struct());
                this.createParam('UseRegExp', false, 'Constraint', biotracs.core.constraint.IsBoolean());
                this.createParam('RemoveNonAnnotated', false, 'Constraint', biotracs.core.constraint.IsBoolean());
                this.createParam('SampleSelectorQuery', '', 'Constraint', biotracs.core.constraint.IsText('IsScalar', false));
                %this.createParam('HasUpdateAnnotations', false, 'Constraint', biotracs.core.constraint.IsBoolean());
          end

	 end
	 
	 % -------------------------------------------------------
	 % Protected methods
	 % -------------------------------------------------------
     
     methods(Access = protected)
  
   
     end

end
