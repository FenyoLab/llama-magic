/*****
**
** functions_ajax.js - a compendium of javascript functions and variables to add some
**	AJAX functionality to the GPM.
**
** Theory of operation: this script will be included in scripts where it is needed. 
**	There is a single ajax_init() function which can be called with various values 
**	to control what is loaded.
**
** Version 2008.06.12 - modified all functions that made AJAX calls to use a global array 
**	to store the calls so that they can be canceled by an external function.  This 
**	is in response to IE not aborting pending AJAX calls when a link on the page was 
**	clicked; it would wait for all AJAX calls to return, and then follow the link.
**
*****/
/*****
**
** GLOBALS
**
*****/

var REQ_OBJECTS = new Array();  // global holder for request objects; key is an int

/*****
**
** FUNCTIONS
**
*****/

function ajax_get_request_object() {
/*
** this function will return a request object of the proper type for the client browser
*/

	var http_request = false;

	if (window.XMLHttpRequest) { // Mozilla, Safari,...

		http_request = new XMLHttpRequest();

		if (http_request.overrideMimeType) {

			http_request.overrideMimeType('text/xml');

		}  // end if

	} else if (window.ActiveXObject) { // IE

		try {  // try to make object

			http_request = new ActiveXObject("Msxml2.XMLHTTP");

		} catch (e) {  // if it fails, try method two

			try {

				http_request = new ActiveXObject("Microsoft.XMLHTTP");

			} catch (e) {}

		}  // end catch

	}  // end else

	if (!http_request) {  // if request object not set

		return false;

	} else { // return object

		return http_request;

	} // end else

} // end ajax_get_request_object()

function ajax_begin_pep_count(id_start, id_end) {
/*
** this script will fire the updates for each <div> tag to be updated.
*/

	var x;
	var div_to_update;
	var seq_value;

	id_end--;  // subtract 1 to make match last id exactly

	for (x=id_start; x <= id_end; x++) {

		div_to_update = 'seqdiv_' + x;
		seq_value = document.getElementById(div_to_update).getAttribute('sequence');
		ajax_count_pep(div_to_update, seq_value);

	}  // end for

}  // end ajax_begin_pep_count()

function ajax_begin_pep_count_batch(id_start, id_end) {
/*
** This function will take a batch of peptides at once, and request that they all be counted.
**	Counting happens in blocks of 20 to get around the URL length limitation of 2048 in 
**	Internet Explorer.  The div ID and accompanying sequence are separated by a = and 
**	the pairs are separated by a |.
**
** Inputs: the start and end numbers associated with the DIVs containing the information.
**
** Outputs: none; calls ajax_count_pep_batch() to do the counting.
*/

	var x;
	var div_name;
	var seq_value;
	var query_data=''; // holder for whole query
	var start_point = id_start;
	var req_ctr=0;
	id_end--; // subtract one for zero indexed div counting
	var blank = '+';

	for (x = id_start; x <= id_end; x++) {

		req_ctr++;
		div_name = 'seqdiv_' + x;
		seq_value = document.getElementById(div_name).getAttribute('sequence');			
		query_data = query_data + div_name + '=' + seq_value + '|';
		blank = blank + div_name + '= |';

		if ((req_ctr == 20) || (x == id_end)) {  // submit update

			query_data = query_data.substring(0, query_data.length - 1);
			ajax_count_pep_batch(query_data,blank);
			query_data = '';  // blank out for next pass
			blank = '+';
			req_ctr = 0;

		}  // end if

	} // end while

}  // end ajax_begin_pep_count_batch()

function ajax_count_pep_batch(delimited_data,blank_data) {
/*
** this function will send the batch request for peptide counts, to try to update more than
**	one at a time.
**
** Inputs: the delimited data; double pipes split data pairs, single pipes split individual pairs
**
** Outputs: the result of the query
*/

	var req = ajax_get_request_object();
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

	var url = '/thegpm-cgi/request_server.pl?target=count_pep&data=' + delimited_data + '&mode=batch';
	req.open('GET', url, true);

	req.onreadystatechange = function() { // define the callback function

		if (req.readyState == 4) {  // got a response, so update with it

	           	if (req.responseText.match("error")) { 

 				ajax_update_innerHTML_batch(blank_data);

           		} else	{

				ajax_update_innerHTML_batch(req.responseText);
			}

		}  // end if

	}  // end anonymous function

	req.send(null);  // send the request

	REQ_OBJECTS.push(ajax_get_request_object());  // add request object to array so it can be aborted if need be

} // end ajax_count_pep_batch()

function ajax_count_pep(div_name, sequence) {
/*
** this function will start the query to count the number of times a peptide has been seen
**
** Inputs: the sequence of the peptide in question.
**
** Output: the result of the query
*/


	var req = ajax_get_request_object();

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

	var url = '/thegpm-cgi/request_server.pl?target=count_pep&data=' + sequence;

	req.open('GET', url, true);

	req.onreadystatechange = function() { // define the callback function

                if (req.readyState == 4) {

	           	if (req.responseText.match("error")) { 
      		        	ajax_update_innerHTML(div_name, ' ');
           		} 
			else	{
				ajax_update_innerHTML(div_name, req.responseText);
			}

                } // end if

        }  // end anonymous function

	req.send(null);  // send the request
	REQ_OBJECTS.push(req);  // add request object to array so it can be aborted if need be

}  // end ajax_count_pep()

function ajax_begin_label_count(id_start, id_end) {
/*
** this function will fire the updates for each <div> tag for the count of label occurrence and 
**	label score.
*/

	var x;
	var div_to_update;
	var label_exp_value;

	id_end--;  // subtract 1 to make match last id exactly

	for (x=id_start; x <= id_end; x++) {

		div_to_update = 'label_' + x;
		label_exp_value = document.getElementById(div_to_update).getAttribute('label_expect');
		ajax_count_label(div_to_update, label_exp_value);

	}  // end for

} // end ajax_begin_label_count()

function ajax_begin_label_count_batch(id_start, id_end) {
/*
** this function will fire the updates for batches of 20 <div> tags for the count of label occurrence
**	and label score.
*/
	var x;
	var div_to_update;
	var label_exp_value;
	var query_data=''; // holder for whole query
	var next_stop = Math.min(id_end, (id_start + 4));  // stopping in 20 or less?
	var update_width = Math.min(4, (id_end - next_stop));  // doing a full 20 or less?
	var start_point = id_start;
	var req_ctr=0;
	id_end--; // subtract one for zero indexed div counting
	var blank = '+';

	for (x = id_start; x <= id_end; x++) {

		req_ctr++;
		div_to_update = 'label_' + x;
		label_exp_value = document.getElementById(div_to_update).getAttribute('label_expect');			
		query_data = query_data + div_to_update + '=' + label_exp_value + '!';
		blank = blank + div_to_update + '= !';

		if ((req_ctr == 4) || (x == id_end)) {  // submit update

			blank = blank + '+';
			query_data = query_data.substring(0, query_data.length - 1);
//			ajax_count_label_batch(query_data,blank);
			ajax_count_label_bucket(query_data, blank);
			query_data = '';  // blank out for next pass
			blank = '+';
			req_ctr = 0;

		}  // end if

	} // end while

}  // end ajax_begin_label_count_batch()

function ajax_count_label_batch(all_data,blank_data) {
/*
** this function requests a batch of counts and expect compares.  The pairs are separated by a 
*	single pipe, and the data parts are separated by a single equal sign.
**
** Inputs: the data as a continuous string.
**
** Outputs: indirect; updates a batch of innerHTML properties.
*/

	var req = ajax_get_request_object();

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

	var url = '/thegpm-cgi/request_server.pl?target=count_label&mode=batch&data=' + all_data;

	req.open('GET', url, true);
	req.onreadystatechange = function() {  // define the callback function

		if (req.readyState == 4) {  // do update

 	           	if (req.responseText.match("error")) { 

				ajax_update_innerHTML_batch(blank_data);

           		} else	{

				ajax_update_innerHTML_batch(req.responseText);

			}  // end else

		}  // end if

	}  // end anonymous function

	req.send(null);
	REQ_OBJECTS.push(req);  // add request object to array so it can be aborted if need be

}  // end ajax_count_label_batch_all()

function ajax_count_label(div_name, label_and_expect) {
/*
** this function requests the count of a specific label
**
** Inputs: the ID of the div tag, and the label and expect value, separated by a double pipe.
**
** Outputs: changes the inner HTML property of the supplied div.
*/

	var req = ajax_get_request_object();

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

	var url = '/thegpm-cgi/request_server.pl?target=count_label&data=' + label_and_expect;

	req.open('GET', url, true);

	req.onreadystatechange = function() { // define the callback function

		if (req.readyState == 4) {

	           	if (req.responseText.match("error")) { 

        		       ajax_update_innerHTML(div_name, ' '); 

           		} else {

				ajax_update_innerHTML(div_name, req.responseText);

			}  // end else

		}  // end if

	}  // end anonymous function

	req.send(null);  // send the request
	REQ_OBJECTS.push(req);  // add request object to array so it can be aborted if need be

} // end ajax_count_label()

function ajax_begin_label_bucket_batch(id_start, id_end) {
/*
** this function will fire the updates for batches of 20 <div> tags for the count of label occurrence
**	and label score.
*/
	var x;
	var div_to_update;
	var label_exp_value;
	var query_data=''; // holder for whole query
	var next_stop = Math.min(id_end, (id_start + 4));  // stopping in 20 or less?
	var update_width = Math.min(4, (id_end - next_stop));  // doing a full 20 or less?
	var start_point = id_start;
	var req_ctr=0;
	id_end--; // subtract one for zero indexed div counting
	var blank = '+';

	for (x = id_start; x <= id_end; x++) {

		req_ctr++;
		div_to_update = 'label_' + x;
		label_exp_value = document.getElementById(div_to_update).getAttribute('label_expect');
	
		query_data = query_data + div_to_update + '=' + label_exp_value + '!';
		blank = blank + div_to_update + '= !';

		if ((req_ctr == 4) || (x == id_end)) {  // submit update

			blank = blank + '+';
			query_data = query_data.substring(0, query_data.length - 1);

			ajax_count_label_bucket(query_data, blank);
			query_data = '';  // blank out for next pass
			blank = '+';
			req_ctr = 0;

		}  // end if

	} // end while

}  // end ajax_begin_label_bucket_batch()

function ajax_count_label_bucket(all_data, blank_data) {
/*
** this function requests a batch of counts and expect compares.  The pairs are separated by a 
*	single pipe, and the data parts are separated by a single equal sign.
**
** Inputs: the data as a continuous string.
**
** Outputs: indirect; updates a batch of innerHTML properties.
*/

	var req = ajax_get_request_object();

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

	var url = '/thegpm-cgi/request_server.pl?target=label_bucket&mode=batch&data=' + all_data;
	req.open('GET', url, true);
	req.onreadystatechange = function() {  // define the callback function

		if (req.readyState == 4) {  // do update

 	           	if (req.responseText.match("error")) { 

				ajax_update_innerHTML_batch(blank_data);


           		} else	{

				ajax_update_innerHTML_batch(req.responseText);

			}  // end else

		}  // end if

	}  // end anonymous function

	req.send(null);
	REQ_OBJECTS.push(req);  // add request object to array so it can be aborted if need be

}  // end ajax_count_label_bucket()

function ajax_update_innerHTML(obj_id, innerHTML) {
/*
** this script will update the inner HTML of a DOM object.
*/

	var object = document.getElementById(obj_id);
	if(!object) {

		return false;

	}  // end if

	object.innerHTML = innerHTML;

}  // end ajax_update_innerHTML()

function ajax_update_innerHTML_batch(delimited_data) {
/*
** this script will update the inner HTML of a batch of DOM objects.  Single pipes delimit the pairs,
**	single pipes delimit the data within each pair.
**
** Inputs: the delimited information
**
** Outputs: indirect; content of the innerHTML property of an object is updated
*/


	if (!delimited_data)	{

		return false;

	}  // end if

	var dirty_data = new Array();
	dirty_data = delimited_data.split('+'); // get middle chunk of data, as it's surrounded by plus signs
	var clean_data = dirty_data[1];

	var pair_array = new Array();  // holds all name/value pairs
	var data_array = new Array();  // holds one split name/value pair
	var object = '';  // DOM object to be updated
	var x = 0;  // counter
	if(!clean_data)	{
		return false;
	}
	if (clean_data.indexOf('!') != -1) {  // check for exclamation point separator for gi|nnn|-type labels
/*
** split on exclamation points instead of a pipe, for input from ajax_count_label_batch
*/

		pair_array = clean_data.split('!');

	} else {  // split normally

		pair_array = clean_data.split('|');  // split on each single pipe

	}  // end else

	for (x=0; x < pair_array.length; x++) {  // for each name/value pair

		data_array = pair_array[x].split('=');  // split into a (name, value) array

		object = document.getElementById(data_array[0]);
		if(object)	{
			object.innerHTML = data_array[1];  // update the content
		}

	}  // end for	

}  // end ajax_update_innerHTML_batch()

function ajax_change_object_vis(obj_name) {
/*
** This function will switch the visibility of the named object.
**
** Inputs: the DOM name of the object to switch.
**
** Outputs: none.
*/

	var obj_to_change = document.getElementById(obj_name);
	var current_vis = obj_to_change.style.display;
	var default_vis = 'block';  // default visibile style

	if (current_vis == 'none') { // make visible

		obj_to_change.style.display = default_vis;

	} else { // make invisible

		obj_to_change.style.display = 'none';

	}  // end else

	var header_img_name = obj_name + '_img';

	if (obj_to_change.style.display == 'none') {  // switch to collapsed pic

		document.getElementById(header_img_name).src = '/pics/ffa_collapsed.gif';

	} else {  // switch to expanded pic

		document.getElementById(header_img_name).src = '/pics/ffa_expanded.gif';

	}  // end else

} // end ajax_change_object_vis()

function ajax_create_cookie(cookie_name, value, expire_days) {
/*
** This function will create a cookie with javascript. the cookies will be used to 
**	store the viewing preferences on some pages of the GPM.
**
** Inputs: the name of the cookie, the value of the cookie, and the number of days 
**	until the cookie exipres.
**
** Outputs: a cookie.
*/

	var date = new Date();
	var expires='';

	if (!expire_days) { // if no day is defined

		expire_days = 30; // set to expire in 30 days

	}  // end if

	date.setTime(date.getTime()+(expire_days*24*60*60*1000));
	expires = "; expires="+date.toGMTString();

	document.cookie = cookie_name + "=" + value + expires + "; path=/; domain=thegpm.org";

} // end ajax_create_cookie()

function ajax_read_cookie(cookie_name) {
/*
** this function will return the data portion of a set cookie.
**
** Inputs: the name of a cookie to read.
**
** Outputs: the data value, or null on failure.
*/

	var return_value = null;  // default return value
	var nametag = cookie_name + '=';
	var cookie_array = document.cookie.split(';');
	var x;

	for(x=0; x < cookie_array.length; x++) {

		var one_cookie = cookie_array[x];
		while (one_cookie.charAt(0)==' ') one_cookie = one_cookie.substring(1,one_cookie.length);

		if (one_cookie.indexOf(nametag) == 0) { // return the value

			return one_cookie.substring(nametag.length, one_cookie.length);

		} // end if

	}  // end for

	return null;  // not found, so return null
	
}  // end ajax_read_cookie()

function ajax_store_display_prefs(page_to_store) {
/*
** This function will store the preferences for a page registered below. It will 
**	set a cookie with some known attribute=value pairs, listed below.
**
** Inputs: the name of the page for which to store preferences.
**
** Outputs: a cookie being set to save the preferences.
*/

	var value_to_store = '';  // value to store as a cookie
	var display_value = 1;  // holder for display value of an element

	if (page_to_store == 'protein') {  // munge together the appropriate content

		var display_value = 1;

		if (document.getElementById('sequence_table').style.display == 'none') {

			display_value = 0;

		}  // end if

		value_to_store = 'sequence_table='+display_value+'&';

		display_value++;

		if (document.getElementById('spectrum_table').style.display == 'none') {

			display_value = 0;

		}  // end if

		value_to_store = value_to_store + 'spectrum_table='+display_value+'&';

		display_value++;

		if (document.getElementById('column_notes').style.display == 'none') {

			display_value = 0;

		}  // end if

		value_to_store = value_to_store + 'column_notes='+display_value+'&';

		display_value++;

		if (document.getElementById('report_section').style.display == 'none') {

			display_value = 0;

		}  // end if

		value_to_store = value_to_store + 'report_section='+display_value+'&';

		if (document.getElementById('n_or_omega_n').checked == true) {

			value_to_store = value_to_store + 'n_or_omega=n';

		} else {

			value_to_store = value_to_store + 'n_or_omega=omega';

		}  // end if

		ajax_create_cookie('protein', value_to_store, 30);

	}  // end if

	ajax_update_innerHTML('display_prefs_msg', "Preferences Saved");  // display status message
	setTimeout("ajax_update_innerHTML('display_prefs_msg', '')",3000); // clear status message

}  // end ajax_store_display_prefs()

function ajax_load_display_prefs(page_to_load) {
/*
** This function will try to load any saved display preferences for a given page
**
** Inputs: the name of the page to check.
**
** Outputs: none direct; will modify the way the page is rendered.
*/

	var cookie_data = '';  // will hold data read from cookie

	if (page_to_load == 'protein') {

		cookie_data = ajax_read_cookie('protein');  // try to read cookie

		if (cookie_data != null) {  // if read successfully

			var data_array = cookie_data.split('&');  // split on &
			var x;  // counter

			for (x=0; x < data_array.length; x++) {  // for each data bit

				var element_vis_arr = data_array[x].split('=');  // split on =
				var element = element_vis_arr[0];  // element name first
				var vis_to_set = element_vis_arr[1];
				var vis;

				if (element == 'n_or_omega') {  // dea with radio button

					var btn_name = 'n_or_omega_' + vis_to_set;
					document.getElementById(btn_name).checked = true;

				} else {

					if (vis_to_set == 0) {  // hide it

						vis = 'none';

					} else {  // show it

						vis = 'block';

					}  // end else
					document.getElementById(element).style.display = vis;

					var pic_name = element + "_img";

					if (vis == 'none') {  // change the pic as well as the text display

						document.getElementById(pic_name).src = '/pics/ffa_collapsed.gif';

					} else {

						document.getElementById(pic_name).src = '/pics/ffa_expanded.gif';

					}  // end if

				}  // end else

			}  // end for

		}  // end if

	}  // end if

}  // end ajax_load_display_prefs()

function ajax_clear_display_prefs(to_clear) {
/*
** This function will clear the stored preferences for the given page.
**
** Inputs: the page to be cleared.
**
** Outputs: a status message in the display_prefs_msg span showing success or failure.
*/

	ajax_create_cookie(to_clear, "", -1);  // clear the cookie
	ajax_update_innerHTML('display_prefs_msg', "Preferences Cleared");
	setTimeout("ajax_update_innerHTML('display_prefs_msg', '')",3000);

}  // end ajax_clear_display_prefs()

function ajax_abort_req_objects() {
/*
** This function will cancel all outstanding requests when a link is clicked; used by IE 
**	browsers for link onclick events, so clicks cancel pending AJAX calls.
**
** Inputs: none.
**
** Outputs: none.
*/
	if(window.XMLHttpRequest == null)	{
		return;
	}
	var obj;
	for (obj in REQ_OBJECTS) {
		if(REQ_OBJECTS[obj] != null && REQ_OBJECTS[obj].abort != null)	{
			REQ_OBJECTS[obj].abort();
		}

	}  // end for

}  // end ajax_abort_req_objects();
