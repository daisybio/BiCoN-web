var currentTab = 0; // Current tab is set to be the first tab (0)
showTab(currentTab); // Display the current tab

function showTab(n) {
    // This function will display the specified tab of the form ...
    var x = document.getElementsByClassName("tab");
    x[n].style.display = "block";
    // ... and fix the Previous/Next buttons:
    if (n === 0) {
        document.getElementById("prevBtn").style.display = "none";
    } else {
        document.getElementById("prevBtn").style.display = "inline";
    }
    if (n === (x.length - 1)) {
        document.getElementById("nextBtn").innerHTML = "Submit";
    } else {
        document.getElementById("nextBtn").innerHTML = "Next";
    }
    // ... and run a function that displays the correct step indicator:
    fixStepIndicator(n)
}

function nextPrev(n) {
    // This function will figure out which tab to display
    var x = document.getElementsByClassName("tab");
    // Exit the function if any field in the current tab is invalid:
    if (n === 1 && !validateForm()) return false;
    // Hide the current tab:
    x[currentTab].style.display = "none";
    // Increase or decrease the current tab by 1:
    currentTab = currentTab + n;
    // if you have reached the end of the form... :
    if (currentTab >= x.length) {
        //...the form gets submitted:
        document.getElementById("regForm").submit();

        // Disable buttons and display waiting page
        document.getElementById("tab-submit").style.display = "block";
        document.getElementById("prevBtn").style.display = "none";
        document.getElementById("nextBtn").style.display = "none";

        return false;
    }
    // Otherwise, display the correct tab:
    showTab(currentTab);
}

function validateForm() {
    // This function deals with validation of the form fields
    var x, y, i, valid = true;
    x = document.getElementsByClassName("tab");
    y = x[currentTab];

    // Check the expression tab
    if (y.classList.contains("tab-expression")) {
        // Check if no option has been selected, set to false
        if ($('input[name=expression-data]:checked').length === 0) {
            valid = false;
            $('#expression_err').addClass('d-block')
            // If upload own is selected, check if file is also selected
        } else if (document.getElementById('expression-data_3').checked
            && document.getElementById("expression-data-filename").value === "") {
            valid = false;
            $('#expression_file_err').addClass('d-block')

        } else {
            $('#expression_err').removeClass('d-block')
            $('#expression_file_err').removeClass('d-block')

        }

        // If predefined option is used, disable metadata selection options
        if (!document.getElementById('expression-data_3').checked) {
            // Mute text
            $('#metadata_options').addClass('text-muted')
            // Disable radios and select yes
            $('input[name=use_metadata]').attr("disabled", true);
            $('input[name=use_metadata]').val(['yes'])
            $('#use_metadata_info').addClass('d-block')
        } else {
            // Unmute text
            $('#metadata_options').removeClass('text-muted')
            // Enable radios and select no
            $('input[name=use_metadata]').attr("disabled", false);
            $('input[name=use_metadata]').prop('checked', false);
            $('#use_metadata_info').removeClass('d-block')
        }

        // Check the ppi tab
    } else if (y.classList.contains("tab-ppi")) {
        // Check if no option has been selected, set to false
        if ($('input[name=ppi-network]:checked').length === 0) {
            valid = false;
            $('#ppi_err').addClass('d-block')
            // If upload own is selected, check if file is also selected
        } else if (document.getElementById('ppi-network_5').checked
            && document.getElementById("ppi-network-filename").value === "") {
            valid = false;
            $('#ppi_file_err').addClass('d-block')

        } else {
            $('#ppi_err').removeClass('d-block')
            $('#ppi_file_err').removeClass('d-block')

        }

        // Metadata tab
    } else if (y.classList.contains("tab-metadata")) {
        // If metadata is not selected or disabled do nothing
        if ($('input[name=use_metadata]:checked', '#regForm').val() === 'yes'
            && !$('input[name=use_metadata]').is('[disabled=disabled]')) {
            // Check if column name is given
            if (!$('#survival-col').val().trim()) {
                $('#survival-col_err').addClass('d-block')
                valid = false
            } else {
                $('#survival-col_err').removeClass('d-block')
            }

            // Check if upload file is given
            if (document.getElementById("survival-metadata").value === "") {
                $('#survival-metadata_err').addClass('d-block')
                valid = false
            } else {
                $('#survival-metadata_err').removeClass('d-block')
            }
        }


        // Check algorithm parameter tab
    } else if (y.classList.contains("tab-parameters")) {
        const L_g_min = $('#L_g_min').val();
        const L_g_max = $('#L_g_max').val();

        // Check if L_g_min is a number
        if (!isNumber(L_g_min)) {
            valid = false;
            $('#L_g_min_err').addClass('d-block')
        } else {
            $('#L_g_min_err').removeClass('d-block')
        }

        // Check if L_g_max is a number
        if (!isNumber(L_g_max)) {
            valid = false;
            $('#L_g_max_err').addClass('d-block')
        } else {
            $('#L_g_max_err').removeClass('d-block')
        }

        // Check if min is smaller or equal to max
        if (valid) {
            if (parseInt(L_g_min) > parseInt(L_g_max)) {
                $('#cluster_param_error').addClass('d-block')
                valid = false;
            } else {
                $('#cluster_param_error').removeClass('d-block')
            }
        }

        // If advanced options are selected
        if ($('input[name=clustering_advanced]:checked', '#regForm').val() === 'yes') {
            const gene_set_size = $('#gene_set_size').val();
            const nbr_iter = $('#nbr_iter').val();

            // Check if gene_set_size is a number
            if (!isNumber(gene_set_size)) {
                valid = false;
                $('#gene_set_size_err').addClass('d-block')
            } else {
                $('#gene_set_size_err').removeClass('d-block')
            }

            // Check if nbr_iter is a number
            if (!isNumber(nbr_iter)) {
                valid = false;
                $('#nbr_iter_err').addClass('d-block')
            } else {
                $('#nbr_iter_err').removeClass('d-block')
            }


        }


    }
    if (valid) {
        document.getElementsByClassName("step")[currentTab].className += " finish";
    }
    return valid; // return the valid status
    // return true; // return the valid status
}

function isNumber(n) {
    return !isNaN(parseFloat(n)) && !isNaN(n - 0)
}


function fixStepIndicator(n) {
    // This function removes the "active" class of all steps...
    var i, x = document.getElementsByClassName("step");
    for (i = 0; i < x.length; i++) {
        x[i].className = x[i].className.replace(" active", "");
    }
    //... and adds the "active" class to the current step:
    x[n].className += " active";
}