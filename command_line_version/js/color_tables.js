var tables = document.getElementsByTagName("table");


for(var i = 0; i < tables.length; i++) {
    var table = tables[i];

    //get all the rows with tag "td"
    var tds = table.getElementsByTagName("td"); 

    for (var j = 0; j < tds.length; j++) {

	var text = tds[j].innerHTML;
	
	//Check what is the last character of a cell, to assign the correct color for matching ions
	var reg = /Mclr$/.test(text)
	if (reg) {

		//Use grey for M ions
		tds[j].style.backgroundColor = "#abb2b9"
		
		var newValue = text.split("_")
		tds[j].innerHTML = newValue[0]
    	}

	var reg = /[a|w]clr$/.test(text)
	if (reg) {

		//Use green for a/w series
		tds[j].style.backgroundColor = "#33CC00"
		
		var newValue = text.split("_")
		tds[j].innerHTML = newValue[0]
    	}

	var reg = /[x|b]clr$/.test(text)
	if (reg) {

		//Use blue for x/b series
		tds[j].style.backgroundColor = "#0099FF"
		
		var newValue = text.split("_")
		tds[j].innerHTML = newValue[0]
    	}

	var reg = /[y|c]clr$/.test(text)
	if (reg) {

		//Use pink for y/c series
		tds[j].style.backgroundColor = "#FF33FF"
		
		var newValue = text.split("_")
		tds[j].innerHTML = newValue[0]
    	}

	var reg = /[z|d]clr$/.test(text)
	if (reg) {

		//Use yellow for z/d series
		tds[j].style.backgroundColor = "#FFCC00"
		
		var newValue = text.split("_")
		tds[j].innerHTML = newValue[0]
    	}

	var reg = /[X]$/.test(text)
	if (reg) {

		//Use cyan for brown for CMC specific series
		tds[j].style.backgroundColor = "#964B00"
		
		var newValue = text.split("_")
		tds[j].innerHTML = newValue[0]
    	}



	var reg = /decoy$/.test(text)
	if (reg) {

		//Use red to highlight decoys in the main table
		tds[j].style.backgroundColor = "#FFC4C4"
		tds[j-14].style.backgroundColor = "#FFC4C4"
		
		
    	}
 

 	}

    
}


