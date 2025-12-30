//////////////////////////// INFORMATION //////////////////////////////
//
//A macro to process root files straight out of LArSoft.
//
//This macro will accept input root files straight from LArSoft and process it to be ready for the next step in the analysis chain.
//
//The output produced will be multiple root files containing all the information for each background category.
//
//Background categories include; NC Deltas, Pi0's, Nues.
//NC deltas will exist from their own simulation so this code will accept 2 root files:
//		1. One containing raw NC Delta truth Information
//		2. One containing raw Pi0 and Nue truth information
//
///////////////////////////////////////////////////////////////////////

///////////////////////// INFORMATION /////////////////////////
//
//This macro sorts through the full root files for the 3 background categories (Pi0, NCDeltas, Nue's), and outputs a root file containing only needed information for the next step 
//
//This macro takes 3 input .root files:
//
//		1. A .root file containing all NC Delta 
//		   MC truth information.
//
//		2. A .root file containing all Pi0
//		   MC truth information.
//
//		3. A .root file containing all Nue
//		   MC truth information
//
//This macro creates 1 output root file: 
//		1. A .root file containing three MC truth
//		   spectra corresponding  to Pi0's,
//		  NC Deltas, and Nue's.
//
/////////////////////////////////////////////////////////////////  

