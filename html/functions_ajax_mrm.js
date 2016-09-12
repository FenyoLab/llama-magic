/*****
**
** functions_ajax_mrm.js - a compendium of javascript functions and variables to add some
**	AJAX functionality to the GPM.
**
** Theory of operation: this script will be included in scripts where it is needed. 
**	There is a single ajax_init() function which can be called with various values 
**	to control what is loaded.
**
** version 2008.06.12 - inaugural version!
** version 2009.04.27 - added functions for the omega calculations with respect to 
**	charge state.
** version 2009.05.27 - changed capital omega formatting to strictly two decimal 
**	places instead of three or four.
**
*****/

/*****
**
** GLOBALS
**
*****/

var SELECTED_PROTEIN = '';  // global to hold the id of the currently selected protein
var PROTEIN_SCORES = new Array();  // global to hold scores by accession number; populated by a function below
var CHARGE_OMEGAS =  new Array();  // global to hold counts by differences in charge
var PEP_COUNT_BY_DIV = new Array();  // peptide totals by div name
var PEP_COUNT_BY_CHARGE = new Array();
var OMEGA_CHARGE_PEP_COUNT_REMAINING;  // global count of peptides to be counted
var OMEGA_SCORES = new Array();  // scores for each element
var PEP_FREQ_BY_CHARGE = new Array();  // two-dimensional array keyed by peptide sequence, then charge.
var PRO_CHARGE_TOTALS = new Array();  // count by charge for each div
var PRO_TOTAL_OBS = 0;
var FIRST_DIV = 0;  // global of first div counter
var LAST_DIV = 0;  // global of last div counter
var Z_TOT = new Array();  // holds totals by charge

/*****
**
** FUNCTIONS
**
*****/

function ajax_get_pro_desc(id, target_id) {
/*
** This script will query and return the description of the protein identified by the given label.
*/

	var req = ajax_get_request_object();
	var div_to_update = document.getElementById(id);
	var label = div_to_update.getAttribute('label');

/*
** in order for links to work properly during pending AJAX requests, IE must abort all transfers when 
**	 a link is clicked.
*/

	if (window.ActiveXObject) {  // do only for IE

		var i = 0;

		for (i=0; i<document.links.length; i++) {  // for each link in this document
			var str = document.links[i].getAttribute("href"); 
			if (str.substring(11,0) != "javascript:" && document.links[i].onclick == null) {  // if not a js link

				document.links[i].onclick = ajax_abort_req_objects;

			} // end if

		} // end for

	} // end if

	var url = '/thegpm-cgi/request_server.pl?target=pep_desc&data=' + label;

	req.open('GET', url, true);

	req.onreadystatechange = function() { // define the callback function

		if (req.readyState == 4) {  // got a response, so update with it

			ajax_set_span_decoration(id, "1px dashed #000000");  // set decoration on 

	           	if (req.responseText.match("error")) {  // update with "(n/a)"

 				ajax_update_innerHTML(id, target_id, "(n/a)");

           		} else	{  // update with what was returned

				ajax_update_innerHTML(id, target_id, req.responseText);

			}  // end else

		}  // end if

	}  // end anonymous function

	req.send(null);  // send the request

	REQ_OBJECTS.push(ajax_get_request_object());  // add request object to array so it can be aborted if need be

}  // end ajax_get_pro_desc()

function ajax_set_span_decoration(obj_id, span_style) {
/*
** This function will set some styling on an object to denote which protein has been clicked
*/

	var obj;

	if (SELECTED_PROTEIN != '') {  // unset style of currently selected protein

		obj = document.getElementById(SELECTED_PROTEIN);
		obj.style.borderBottom = 'none';

	}  // end if

	SELECTED_PROTEIN = obj_id;
	obj = document.getElementById(obj_id);
	obj.style.borderBottom = span_style;

}  // end ajax_set_span_decoration()

function ajax_update_protein_details(id, target_id) {
/*
** This function will update the content of the "details" box that pops up when a protein is clicked 
**	in the peak_details.pl window.
**
** NOTE: this function deprecated, and may be removed
*/

	var obj = document.getElementById(target_id);
	obj.style.position = 'relative';

}  // end ajax_update_protein_details()

function ajax_begin_omega_score(start, end) {
/*
** gets all the label/sequence data for a given set of omega score span tags
*/

	var x;
	var obj;
	var span_name;
	var label;
	var seq;

	for (x=start; x<end; x++) {  // get list of label/sequence tuples

		span_name = 'omega_' + x;
		obj = document.getElementById(span_name);
		label = obj.getAttribute('label');
		seq = obj.getAttribute('sequence');
		ajax_get_omega_scores(span_name, seq);  // get the score for this one peptide

	}  // end for

}  // end ajax_begin_omega_score()

function ajax_get_omega_totals(label) {
/*
** gets the peptide counts per charge state for a specific protein
**
** PLEASE NOTE: this function has been deprecated, and will be removed.
*/

	var req = ajax_get_request_object();
	var url = '/thegpm-cgi/request_server_mrm.pl?target=pro_charge_totals&data=' + label;

	req.open('GET', url, false);  // needs to be synchronous to ensure the omega calculations are right
	req.send(null);  // send the request
	ajax_set_protein_charge_counts(req.responseText);

}  // end ajax_get_omega_totals()

function ajax_calculate_big_omega(start, end) {
/*
** This function calculate the big omega proportions for the charge states 1-3 in the 
**	peptides displayed on protein.pl.
*/

	var x;
	var div;
	var seq_z_score = new Array();  // keyed as seq_z_score[seq][z]
	var big_omegas = new Array();
	var big_omega_div_name = 'protein_omega_' + SELECTED_PROTEIN;
	var big_omega_str = '';
	var one_omega;
	var one_z;
	var one_seq;
	var rounded_score;

	for (x=start; x<end; x++) {  // for each div with information

		div = document.getElementById('seqdiv_' + x);
		one_z = div.getAttribute('charge');
		one_seq = div.getAttribute('sequence');

		if (one_z < 4) {  // try to record

			one_omega = div.innerHTML;

			if (!isNaN(one_omega)) {  // record

				if (typeof(seq_z_score[one_seq]) == 'undefined') {  // record

					seq_z_score[one_seq] = new Array();
					seq_z_score[one_seq][one_z] = parseFloat(one_omega);

				} else {

					if (typeof(seq_z_score[one_seq][one_z]) == 'undefined') {  // record

						seq_z_score[one_seq][one_z] = parseFloat(one_omega);

					}  // end if

				}  // end else

			}  // end if

		}  // end if

	}  // end for

	for (one_seq in seq_z_score) {  // for each seq

		for (one_z in seq_z_score[one_seq]) {  // for each charge

			one_omega = seq_z_score[one_seq][one_z];

			if (!isNaN(one_omega)) {  // add

				if (typeof(big_omegas[one_z]) == 'undefined') {  // define

					big_omegas[one_z] = parseFloat(one_omega);

				} else {  // increment

					big_omegas[one_z] += parseFloat(one_omega);

				}  // end else

			}  // end if

		}  // end for

	}  // end for

	for (x=1; x<4; x++) {

		rounded_score = Math.round(big_omegas[x] * 100)/100;
		if(isNaN(rounded_score))	{
			rounded_score = 0;
		}
		big_omegas[x] = rounded_score;

	}  // end for
	big_omega_str = big_omega_str + '&Omega; = ';
	big_omega_str = big_omega_str + big_omegas[1] + '<sup>1</sup>|';
	big_omega_str = big_omega_str + big_omegas[2] + '<sup>2</sup>|';
	big_omega_str = big_omega_str + big_omegas[3] + '<sup>3</sup>';

//	big_omega_str = big_omega_str.substring(0, big_omega_str.length - 6);  // remove trailing space
	ajax_update_innerHTML(big_omega_div_name, big_omega_str);  // update the div

}  // end ajax_calculate_big_omega()

function ajax_set_protein_charge_counts(response) {
/*
** records the peptide counts per charge state into the global variable PRO_CHARGE_TOTALS
*/

	var resp_arr = response.split('+');
	var count_arr = resp_arr[1];  // data in the middle chunk
	var z_total_arr = count_arr.split('^');
	var big_omega = new Array();  // top of the column display
	var big_omega_text = '&Omega;&nbsp;';
	var div_name = 'protein_omega_' + SELECTED_PROTEIN;  // the div to hold the overall frequency distribution
	var x;
	var y;
	var value;

	for (x=0; x<z_total_arr.length; x++) {  // for each element

		y = x+1;
		PRO_CHARGE_TOTALS[y] = z_total_arr[x];
		PRO_TOTAL_OBS += parseInt(z_total_arr[x]);

	}  // end for

	for (x=0; x < z_total_arr.length; x++) {  // for each element

		y = x+1;
		value = PRO_CHARGE_TOTALS[y] / PRO_TOTAL_OBS;
		big_omega[y] = Math.round(value * 100) / 100;

	}  // end for

	for (x=0; x<z_total_arr.length; x++) {  // for each element

		y = x+1;
		big_omega_text = big_omega_text + 'z<sup>' + y + '</sup>=' + big_omega[y] + '&nbsp;';

	}  // end for

	big_omega_text = big_omega_text.substring(0, big_omega_text.length - 6);  // remove trailing &nbsp;
	ajax_update_innerHTML(div_name, big_omega_text);

}  // end ajax_set_protein_charge_counts()

function ajax_begin_omega_charge_score_batch(start, end) {
/*
** gets all the label/sequence/charge data for a given set of omega score span tags
*/

	OMEGA_CHARGE_PEP_COUNT_REMAINING = end - start;  // the number of peptides to still be returned
	FIRST_DIV = start;
	LAST_DIV = end;
	var x;
	var obj;
	var span_name;
	var label;
	var seq;
	var charge;
	var data_string = '';
	var ctr=1;
	var res = '';
	var res_array;
	var data = '';
	var omega_totals;
	var y = 1;
	span_name = 'seqdiv_' + start;
	obj = document.getElementById(span_name);
	label = obj.getAttribute('label');
	SELECTED_PROTEIN = label;  // set global protein name

	for (x=start; x<end; x++) {  // get the list of label/seq/charge tuples

		span_name = 'seqdiv_' + x;
		obj = document.getElementById(span_name);

		seq = obj.getAttribute('sequence');
		charge = obj.getAttribute('charge');

		data = data + span_name + '=' + seq + '^' + charge + '!';

		if (((y % 5) == 0) || (y == (end - start))) {  // submit batch

			data = data.substring(0, data.length - 1);  // remove trailing separator
			ajax_get_pep_omega_score(data, SELECTED_PROTEIN);
			data = '';

		}  // end if

		y++;

	}  // end for

}  // end ajax_begin_omega_charge_score_batch()

function ajax_get_pep_omega_score(tuples, label) {
/*
** gets the omega score for a peptide in the requested protein
*/

	var url = '/thegpm-cgi/request_server_mrm.pl?target=pep_omega_charge&data=' + label + '^^' + tuples;
	var req = ajax_get_request_object();
	req.open('GET', url, true);

	req.onreadystatechange = function() { // define the callback function

		if (req.readyState == 4) {  // got a response, so update with it

			res = req.responseText;
			res_array = res.split('+');
			data = res_array[1];  // get middle chunk
			ajax_set_pep_omega_score(data);  // set returned batch of data

		}  // end if

	}  // end anonymous function

	req.send(null);  // send the request

	REQ_OBJECTS.push(ajax_get_request_object());  // add request object to array so it can be aborted if need be

}  // end ajax_get_pep_omega_score()

function ajax_set_pep_omega_score(data) {
/*
** sets the omega scores for specific peptides
*/

	var tuples = data.split('!');
	var x;
	var y;
	var div_data_arr;
	var count_z_arr;
	var omega_score;
	var div;
	var z;
	var n;
	var pro_charge_totals_str = '';
	var pro_charge_totals_arr;

	for (x in tuples) {  // for each tuple returned

		div_data_arr = tuples[x].split('=');
		PRO_CHARGE_TOTALS[div_data_arr[0]] = div_data_arr[1];  // store count/charge by div
		OMEGA_CHARGE_PEP_COUNT_REMAINING--;  // decrement global count of peptides left to update

	}  // end foreach

	if (OMEGA_CHARGE_PEP_COUNT_REMAINING == 0) {  // update the display

		pro_charge_totals_str = ajax_get_pro_charge_totals(SELECTED_PROTEIN);
		pro_charge_totals_arr = pro_charge_totals_str.split('^');

		for (x=0;x<3;x++) {  // record the totals per charge state for this protein

			y = x+1;	
			Z_TOT[y] = pro_charge_totals_arr[x];	

		}  // end for

		for (x=FIRST_DIV; x<LAST_DIV; x++) {  // update the peptide omegas

			div = 'seqdiv_' + x;
			count_z_arr = PRO_CHARGE_TOTALS[div].split('^');
			omega_score = ajax_calculate_omega(count_z_arr[0], count_z_arr[1]);
			ajax_update_innerHTML(div, omega_score);  // display the omega score

		}  // end for

		ajax_calculate_big_omega(FIRST_DIV, LAST_DIV);  // fill in the big omega score

	}  // end if

}  // end ajax_set_pep_omega_score()

function ajax_calculate_omega(count, charge) {
/*
** calculates the omega score for a given peptide count and charge state
*/

	var score;
	var formatted_score = '';
	var rounded_score;

	if ((count == 0) || (charge > 3) || (isNaN(count))) {  // catch zeros

		formatted_score = '&mdash;';

	} else {

		score = count / Z_TOT[charge];

	}  // end else

	rounded_score = Math.round(score * 1000)/1000;

	if (rounded_score < 0.01) {

		rounded_score = Math.round(score * 10000)/10000;

	}  // end if

	if (rounded_score < 0.001) {

		rounded_score = Math.round(score * 100000)/100000;

	}  // end if

	if(isNaN(rounded_score)) {

		rounded_score = '&mdash;';

	}  // end if

	if (formatted_score != '') {  // a score to report

		return formatted_score;

	} else {

		return rounded_score;

	}  // end else

}  // end ajax_calculate_omega()

function ajax_get_pro_charge_totals(label) {
/*
** gets the observation totals by charge state for a specific protein
*/

	var req = ajax_get_request_object();
	var url = '/thegpm-cgi/request_server_mrm.pl?target=pro_charge_totals&data=' + label;
	var resp;  // response;
	req.open('GET', url, false);  // must be synchronous
	req.send(null);
	resp = req.responseText.split('+');
	return resp[1];  // the middle chunk is what's needed

}  // end ajax_get_pro_charge_totals()

function ajax_get_omega_charge_batch(data) {
/*
** This function will request the total peptide count for a protein and the specific 
**	identification counts for each peptide, taking charge state into account
*/

	var req = ajax_get_request_object();
	var url = '/thegpm-cgi/request_server_mrm.pl?target=pro_pep_charge_omega&data=' + data;
	var res = '';

	req.open('GET', url, true);

	req.onreadystatechange = function() { // define the callback function

		if (req.readyState == 4) {  // got a response, so update with it

			res = req.responseText;
			res_array = res.split('+');
			data = res_array[1];  // get middle chunk
			ajax_record_omega_charge_scores(data);  // set returned batch of data

		}  // end if

	}  // end anonymous function

	req.send(null);  // send the request

	REQ_OBJECTS.push(ajax_get_request_object());  // add request object to array so it can be aborted if need be

}  // end ajax_get_omega_charge_batch()

function ajax_record_omega_charge_scores(data) {
/*
** this function parses the response from the AJAX server and adds the scores to the global
**	variables.
*/

	var data_arr = data.split('!');
	OMEGA_CHARGE_PEP_COUNT_REMAINING -= data_arr.length;
	var div_count_arr;
	var pair;
	var z;
	var obj;
	var x;

	for (x=0; x < data_arr.length;x++) {

		pair = data_arr[x];
		div_count_arr = pair.split('=');

		if (typeof(PEP_COUNT_BY_DIV[div_count_arr[0]]) != "undefined") {

			PEP_COUNT_BY_DIV[div_count_arr[0]] += parseInt(div_count_arr[1]);

		} else {

			PEP_COUNT_BY_DIV[div_count_arr[0]] = parseInt(div_count_arr[1]);

		}  // end else

		obj = document.getElementById(div_count_arr[0]);
		z = obj.getAttribute('charge');

		if (typeof(PEP_COUNT_BY_CHARGE[z]) != "undefined") {

			PEP_COUNT_BY_CHARGE[z] += parseInt(div_count_arr[1]);

		} else {

			PEP_COUNT_BY_CHARGE[z] = parseInt(div_count_arr[1]);

		}  // end else

	}  // end for

	if (OMEGA_CHARGE_PEP_COUNT_REMAINING == 0) {  // do the calculations for each charge state and display

		ajax_calculate_omega_charge();

	}  // end if

}  // end ajax_record_omega_charge_scores()

function ajax_calculate_omega_charge() {
/*
** uses the global arrays defined at the top and populated by the omega_charge-related functions above 
**	to display the omega scores on a per-charge basis for each peptide
*/

	var div;
	var obj;
	var z;
	var charge_total;
	var one_omega;
	var rounded_score=1;
	var BIG_OMEGA = new Array();
	var OMEGA_TEXT = '';
	var OBS_TOTAL = 0;

	for (div in PEP_COUNT_BY_DIV) {  // for each div

		obj = document.getElementById(div);
		z = obj.getAttribute('charge');
		charge_total = PEP_COUNT_BY_CHARGE[z];

		if (charge_total != 0) {  // catch divide by zero error

			one_omega = PEP_COUNT_BY_DIV[div] / charge_total;

		} else {

			one_omega = 0;

		}  // end else

		rounded_score = Math.round(one_omega * 1000)/1000;

		if (rounded_score < 0.01) {

			rounded_score = Math.round(one_omega * 10000)/10000;

		}  // end if

		if (rounded_score < 0.001) {

			rounded_score = Math.round(one_omega * 100000)/100000;

		}  // end if

		if(isNaN(rounded_score)) {

			rounded_score = '&mdash;';

		}  // end if

		if (z != 4) {  // update normally

			ajax_update_innerHTML(div, '&omega;<sup>' + z + '</sup>=' + rounded_score);  // set one score

		} else {

			ajax_update_innerHTML(div, '&omega;<sup>' + z + '</sup>=&mdash;');  // set one score

		}  // end else

	}  // end for
	
	for (z in PEP_COUNT_BY_CHARGE) {

		OBS_TOTAL += parseInt(PEP_COUNT_BY_CHARGE[z]);

	}  // end for

	for (z in PEP_COUNT_BY_CHARGE) {

		BIG_OMEGA[z] = PEP_COUNT_BY_CHARGE[z] / OBS_TOTAL;
		rounded_score = Math.round(BIG_OMEGA[z] * 1000) / 1000;

		if (rounded_score < 0.01) {

			rounded_score = Math.round(BIG_OMEGA[z] * 10000)/10000;

		}  // end if

		if (rounded_score < 0.001) {

			rounded_score = Math.round(BIG_OMEGA[z] * 100000)/100000;

		}  // end if

		if(isNaN(rounded_score)) {

			rounded_score = '&mdash;';

		}  // end if

		OMEGA_TEXT += '&Omega;<sup>' + z + '</sup>=' + rounded_score + ' ';

	}  // end for

	ajax_update_innerHTML('protein_omega_' + SELECTED_PROTEIN, OMEGA_TEXT);

}  // end ajax_calculate_omega_charge()

function ajax_get_omega_score(id, seq_list) {
/*
** This function will request the total peptide count for a protein and the 
**	specific identification counts for each requested sequence.
*/

	var req = ajax_get_request_object();
	var div_to_update = document.getElementById(id);
	var label = div_to_update.getAttribute('label');
	var url = '/thegpm-cgi/request_server_mrm.pl?target=pro_pep_omegas&data=' + label + '^^' + seq_list;
	req.open('GET', url, true);

	req.onreadystatechange = function() { // define the callback function

		if (req.readyState == 4) {  // got a response, so update with it

	           	if (req.responseText.match("error") || req.responseText.match('^^0^')) {  // update with "(n/a)"

 				ajax_update_innerHTML(id, "(n/a)");

           		} else	{  // update with what was returned

				var score = ajax_compute_omega_score(req.responseText);
				var rounded_score = Math.round(score * 1000)/1000;
				if(score < 0.01)	{
					rounded_score = Math.round(score * 10000)/10000;
				}
				if(score < 0.001)	{
					rounded_score = Math.round(score * 100000)/100000;
				}
				if(isNaN(rounded_score))	{
					rounded_score = 0;
				}
				ajax_update_innerHTML(id, '&omega;=' + rounded_score);  // set one score

			}  // end else

		}  // end if

	}  // end anonymous function

	req.send(null);  // send the request

	REQ_OBJECTS.push(ajax_get_request_object());  // add request object to array so it can be aborted if need be

}  // end ajax_get_omega_score()

function ajax_get_omega_score_batch(data_list) {
/*
** function to request a batch of omega scores.
**
** Input: a delimited list of the form "label"^^(div_id=sequence)!+", where multiple 
**	id=sequence tuples are separated by a !.
**
** Output: calls ajax_update_innerHTML_batch_omega and ajax_set_protein_omega_score; 
**		populates the global array PROTEIN_SCORES[][] with values.
*/


	var data_array = data_list.split('^^');
	var label = data_array[0];
	var tuple_list = data_array[1];
	var tuples = tuple_list.split('!');
	var tuple_array;
	var div_id;
	var one_seq;
	var one_tuple;
	var x;
	var submit_data = '';
	var req = ajax_get_request_object();

	for (x=0; x<tuples.length; x++) {

		one_tuple = tuples[x];
		tuple_array = one_tuple.split('=');
		div_id = tuple_array[0];
		one_seq = tuple_array[1];

		if (submit_data == '') {

			submit_data = label + '^^' + div_id + '=' + one_seq;

		} else {  // concatenate

			submit_data = submit_data + '!' + div_id + '=' + one_seq;

		}  // end else

	}  // end for

	var params = 'target=pro_pep_omegas&mode=batch&data=' + submit_data;
	var url = '/thegpm-cgi/request_server_mrm.pl';

	req.open('POST', url, true);
	req.setRequestHeader("Content-type", "application/x-www-form-urlencoded");
	req.setRequestHeader("Content-length", params.length);
	req.setRequestHeader("Connection", "close");

	req.onreadystatechange = function() { // define the callback function

		if (req.readyState == 4) {  // got a response, so update with it

			var response_array = req.responseText.split('+');
			var score_text = response_array[1];
			var score_array = score_text.split('^^');
			var total = score_array[0];
			var tuple_list = score_array[1];
			var tuple_array;
			var one_id;
			var one_score;
			var one_count;
			var one_tuple;
			var id_count_arr;
			tuple_array = tuple_list.split('!');
			var omega_list = '';
			var rounded_score;
			var x;

			for (x=0; x < tuple_array.length; x++) {

				one_tuple = tuple_array[x];

				if (one_tuple == '') {continue;}  // skip blanks

				id_count_arr = one_tuple.split('=');
				one_id = id_count_arr[0];
				one_count = id_count_arr[1];

				if (typeof(one_count) != 'NaN') {  // do calculation

					one_score = one_count / total;
					rounded_score = Math.round(one_score * 1000)/1000;  // three places
					if(one_score < 0.01)	{
						rounded_score = Math.round(one_score * 10000)/10000;
					}
					if(one_score < 0.001)	{
						rounded_score = Math.round(one_score * 100000)/100000;
					}

				} else {  // set to n/a

					rounded_score = 'n/a';

				}  // end else

/*
** add score to PROTEIN_SCORES double array
*/

				var pro_label = document.getElementById(one_id).getAttribute('label');
				var pep_seq = document.getElementById(one_id).getAttribute('sequence');

				if (PROTEIN_SCORES[pro_label][pep_seq] == 0) {  // set the score

					PROTEIN_SCORES[pro_label][pep_seq] = rounded_score;

				}  // end if

				if (omega_list == '') {

					omega_list = one_id + '=' + rounded_score;

				} else {  // concatenate

					omega_list = omega_list + '!' + one_id + '=' + rounded_score;

				}  // end else

			}  // end for

			ajax_update_innerHTML_batch_omega(omega_list);
			ajax_set_protein_omega_score(pro_label);  // each call deals with the same protein

		}  // end if

	}  // end anonymous function

	req.send(params);  // send the request
	REQ_OBJECTS.push(ajax_get_request_object());  // add request object to array so it can be aborted if need be

}  // end ajax_get_omega_score_batch()

function ajax_begin_omega_score_batch(start_id, end_id) {
/*
** function to batch submit omega score requests on a per-protein basis.
**
** Inputs: the start and end id numbers for the span tags to query for data and update.
**
** outputs: populates the global array PROTEIN_SCORES with keys; calls 
**	ajax_get_omega_scores_batch for each label to populate.
*/

	var x;
	var div_name;
	var seq_value;
	var acc;
	var start_point = start_id;
	var req_ctr = 0;
	end_id--;  // subtract one to make it match a zero-indexed count
	var blank = '+';
	var obj;
	var label_array = new Array();  // associative array of labels

	for (x = start_id; x<= end_id; x++) {

		req_ctr++;
		div_name = 'seqdiv_' + x;
		obj = document.getElementById(div_name);
		acc = obj.getAttribute('label');

		if (typeof(PROTEIN_SCORES[acc]) == 'undefined') {  // make an entry for this protein

			PROTEIN_SCORES[acc] = new Array();

		}  // end if

		seq_value = obj.getAttribute('sequence');

		if (typeof(PROTEIN_SCORES[acc][seq_value]) == 'undefined') {  // make an entry for this peptide

			PROTEIN_SCORES[acc][seq_value] = 0;

		}  // end if

		if (typeof(label_array[acc]) == 'undefined') {  // label not yet seen

			label_array[acc] = acc + '^^' + div_name + '=' + seq_value;

		} else {  // concatenate new sequence in

			label_array[acc] = label_array[acc] + '!' + div_name + '=' + seq_value;

		}  // end else

	}  // end for

	for (label in label_array) {  // for each label

		ajax_get_omega_score_batch(label_array[label]);

	}  // end foreach

}  // end ajax_begin_omega_score_batch()

function ajax_compute_omega_score(response_text) {
/*
** computes the score for a particular protein, based off the total number of peptide 
**	identifications and the share of specific peptides within a specific set of 
**	observations.
**
** Inputs: the response text from the AJAX server in the form "+label^^total^(id=score^)+"
**	where there may be multiple id=score tuples separated from each other by a carat.
**
** Output: the overall omega score for a particular peptide.
*/

	var label;
	var scores;
	var id_total;
	var seq_counts;
	var one_seq;
	var one_count;
	var score = 0;  // the overall score for this protein
	var x;
	var arr;

	arr = response_text.split('+');  // we want the second chunk
	var score_text = arr[1];

	arr = score_text.split('^^');
	label = arr[0];
	scores = arr[1];
	arr = scores.split('^');
	id_total = arr[0];  // total identifications for this label isolated
	var arr2;

	for (x=1; x<arr.length; x++) {  // for each sequence's total

		arr2 = arr[x].split('=');
		one_seq = arr2[0];
		one_count = arr2[1];
		score += (one_count/id_total);

	}  // end for

	if (score > 0) {

		return score;

	} else {  // error catch

		return 'n/a';

	}  // end else

}  // end ajax_compute_omega_score()

function ajax_set_protein_omega_score(label) {
/*
** function to get the total score for a given protein
**
** Inputs: the label for which to retrieve a score.
**
** Outputs: calls ajax_update_innerHTML.
*/

	var sum = 0;
	var x;
	var protein_omega_id = 'protein_omega_' + label;
	var one_score;
	var score_array = new Array();
	score_array = PROTEIN_SCORES[label];
	var seq;
	var rounded_sum;

	for (seq in score_array) {

		one_score = score_array[seq];

		if (typeof(one_score) == 'number') {

			sum += one_score;

		}  // end if

	}  // end for

	rounded_sum = Math.round((sum * 100))/100;
	if(sum < 0.1)	{
		rounded_sum = Math.round(sum * 1000)/1000;
	}
	if(sum < 0.01)	{
		rounded_sum = Math.round(sum * 10000)/10000;
	}
	if(sum < 0.001)	{
		rounded_sum = Math.round(sum * 100000)/100000;
	}
	if(isNaN(rounded_sum))	{
		rounded_sum = 0;
	}
	ajax_update_innerHTML(protein_omega_id, rounded_sum);

}  // end ajax_set_protein_omega_score()

function ajax_update_innerHTML_batch_omega(data_pairs) {
/*
** function to update several ids at a time.
**
** Inputs: a list of id=score pairs, delimited by !.
**
** Outputs: calls ajax_update_innerHTML().
*/

	var tuples = data_pairs.split('!');
	var one_pair;
	var id;
	var value;
	var pair_array;
	var x;

	for (x=0; x<tuples.length; x++) {

		one_pair = tuples[x];
		pair_array = one_pair.split('=');
		id = pair_array[0];
		value = pair_array[1];
		if(isNaN(value))	{
			value = 0;
		}
		ajax_update_innerHTML(id, value);  // defined in functions_ajax.js	

	}  // end foreach

}  // end ajax_update_innerHTML_batch_omega()
