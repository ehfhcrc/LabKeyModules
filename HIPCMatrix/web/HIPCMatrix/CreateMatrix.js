/* CreateMatrix.js */

function createMatrixClick(dataRegion)
{
    var schemaName = dataRegion.schemaName;
    var queryName = dataRegion.queryName;
    var selectionKey = dataRegion.selectionKey;
    var cohorts = dataRegion.getColumn("Cohort")
    
    function onlyUnique(value, index, self) { 
    	return self.indexOf(value) === index;
	}

	var uniqueCohorts = cohorts.filter( onlyUnique )

	/* Ensure only one cohort, otherwise names are concatenated and are not found
	by various schema since they are not in the arm_2_cohort table */
	if( uniqueCohorts > 1 ) {
		alert("Please redo selection to use only one cohort.")
	} else {
		window.location = LABKEY.ActionURL.buildURL("HIPCMatrix", "CreateMatrix.view", null, {
	        schemaName: schemaName,
	        queryName: queryName,
	        selectionKey: selectionKey,
	        returnUrl: window.location
    	});
	}
}

    
