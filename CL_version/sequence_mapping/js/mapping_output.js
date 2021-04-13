var tables = document.getElementsByTagName("table");

for(var i = 0; i < tables.length; i++) {
    var table = tables[i];
    var maxSp = 0;

    //get all the rows with tag "td"
    var tds = table.getElementsByTagName("td");

    //get all the rows with tag "th"
    var ths = table.getElementsByTagName("th");

    //Leave only sequence numbers multiple of 10 on the first line
    for (var j = 0; j < ths.length; j++) {

	var value = ths[j].innerHTML;
	
	var reg = /^\d+$/.test(value);
	if (reg) {		
		
		if (value % 10 != 0 && value != 1) {
		
			ths[j].innerHTML = "";	
		} else {
			ths[j].style.fontSize = "8px";
		}
        }

    }
    //Find the max value of the Sp score in the output [maxSp]
    for (var j = 0; j < tds.length; j++) {

	var text = tds[j].innerHTML;
	
	//Find the maximum Sp value for matching occurrences	
	var reg = /^[light|heavy]/.test(text)
	if (reg) {	
		var textSplit = text.split("_");
		
		var Sp = Number(textSplit[2]);        	

       		if (Sp > maxSp) {
			maxSp = Sp;
    		}
	}

    }

    //Modify cells by adding residues names, terminal residues and color coding based on Sp score
    for (var j = 0; j < tds.length; j++) {

	var text = tds[j].innerHTML;
	
	//Check if the cell string starts with 'light' or 'heavy'
	var reg = /^[light|heavy]/.test(text)
	if (reg) {
		
		// Mark the terminal residues for fragment by bolding residue names
		var at = /[\@]/.test(text);
		if (at) {
			tds[j].style.fontWeight = "bold";				
			
			var newStr = text.replace("@", "");
			tds[j].innerHTML = newStr;
		}

		var text = tds[j].innerHTML;
		var ScoreSplit = text.split("_");

		//Mark the cell with the matched residue name	
		tds[j].innerHTML = ScoreSplit[6];

		//Makes the text bigger and the cell grey if a modified base is found
		if (ScoreSplit[6] != 'C' && ScoreSplit[6] != 'G' && ScoreSplit[6] != 'A' && ScoreSplit[6] != 'U' && ScoreSplit[6] != 'X') {

			tds[j].style.fontSize = "10px";
			tds[j].style.border = "1px solid black";			

		} else {
			
			// Unmodified bases are rendered in grey	
			tds[j].style.fontSize = "9px";
			tds[j].style.color = '#3f3f3f';
			tds[j].innerHTML = "";
		}
		
		var Sp = Number(ScoreSplit[2]);
		var ratioSp = Number(Sp / maxSp);
		
		//Color the cell based on the score, with 3 color levels from lower (yellow) to higher (red)
		if (ratioSp <= 0.33) {

			tds[j].style.backgroundColor = "#ffff00";

		} else if (ratioSp > 0.33 && ratioSp < 0.66) {

			tds[j].style.backgroundColor = "#ff9900";

		} else if (ratioSp >= 0.66) {
			tds[j].style.backgroundColor = "#ff0000";
		}
			
    	}
 

    }

}

