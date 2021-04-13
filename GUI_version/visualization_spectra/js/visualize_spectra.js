var title = document.title;
var pics = document.getElementsByTagName("img");

for(var i = 0; i < pics.length; i++) {

	var pic = pics[i];

	//Determine the name of the spectra in the scored_spectra directory
	var spec_id = pic.id;
	var spec_split = spec_id.split("_")
	
	//Add the src to the spectrum with the correct directory in order to visualize it
	pic.src = "scored_spectra_" + title + "/" + spec_split[0] + "_" + spec_split[1] + "_" + spec_split[2] + "_" + spec_split[3] + ".png";


} 
