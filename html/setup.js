
$(document).ready(function()
{
    $( "#dialog" ).dialog( {autoOpen: false, modal:true, closeOnEscape:false, beforeClose: function() { return false; }, title: "Uploading...", draggable:false, resizable:false } );
    
    $('#db_upload_button').button().click(function(e) { $('#dialog').dialog("open") })
    
    $('#ms_upload_button').button().click(function(e) { $('#dialog').dialog("open") })
});

function ec(tId, clickIcon)
{
        dstyle = document.getElementById(tId).style.display;
        if (dstyle == "none")
        {
                document.getElementById(tId).style.display = "";
                document.getElementById(clickIcon).src = "/llama-magic-html/minus.gif";
                document.getElementById(clickIcon).alt = "-";
        }
        else
        {
                document.getElementById(tId).style.display = "none";
                document.getElementById(clickIcon).src = "/llama-magic-html/plus.gif";
                document.getElementById(clickIcon).alt = "+";
        }
}

function switch_right(id, db_id, db_name)
{
    document.getElementById(id).style.display="";
    
    //if ms upload page: fill out text input field and hidden field with db name and id
    if(id == "ms_form")
    {
        document.getElementById("ms_db_name").value=db_name;
        document.getElementById("ms_db_id").value=db_id;
        document.getElementById("ms_db_name").disabled=true;
    }
}

function close_right(id)
{
    document.getElementById(id).style.display="none";
}