let error_counter = 0;
const error_threshold = 100;
let task_status;
let progress_step;
let progress_percent;
let error_message;

function update_status(csrf_token, task_id, poll_url) {
    // Get task_status via AJAX
    $.ajax({
        url: poll_url,
        type: 'POST',
        data: {
            task_id: task_id,
            csrfmiddlewaretoken: csrf_token,
        },
        success: function (result) {
            task_status = result.task_status;
            progress_step = result.progress_step;
            progress_percent = result.progress_percent;
            error_message = result.error_message;
        },

        error: function (error) {
            error_counter++;
        }
    });

    switch (task_status) {
        case 'RUNNING':
            show_only_alert('alert_running');
            $('#progress-bar').html(progress_percent + '% ');
            $('#progress-bar').css('width', progress_percent + '%');

            display_step(progress_step, false);
            break;

        case 'SUCCESS':
            show_only_alert('alert_success');
            $('#progress-bar').html(100 + '% ');
            $('#progress-bar').css('width', 100 + '%');

            display_step('visualize_data', false);
            break;

        // Catched exception
        case 'ERROR':
            show_only_alert('alert_failure');
            $('#progress-bar').html(progress_percent + '% ');
            $('#progress-bar').css('width', progress_percent + '%');

            display_step(progress_step, true);

            $('#error-div').removeClass('d-none').addClass('d-block');
            $('#error-message').html(error_message);
            break;

        // FAILURE means uncatched exception...
        case 'FAILURE':
            show_only_alert('alert_failure');

            style_steps('submitted', false)

            $('#error-div').removeClass('d-none').addClass('d-block');
            $('#error-message').html('An unknown error has occured.\n' +
                'Please contact the administrator of this site and provide the ' +
                'settings you used in this run.\n' +
                'Please also provide this task-id: ' + task_id);

            break;

        case 'PENDING':
            show_only_alert('alert_pending');
            break;
    }


    // Update graph and text while the process is running
    if (!(task_status === 'FAILURE' || task_status === 'SUCCESS' || task_status === 'ERROR'
        || error_counter > error_threshold)) {
        setTimeout(update_status, 2000, csrf_token, task_id, poll_url);
        console.log("Update status")
    }
}

// Color the steps and set right symbols
function style_steps(step_number, failure) {
    jQuery('.process-step ').each(function (index, currentElement) {

        // Finished steps
        if (index < step_number) {
            $(currentElement).removeClass('text-muted').addClass('text-success');
            $(currentElement).removeClass('empty-circle').addClass('check-circle');

            // Current step
        } else if (index === step_number) {
            $(currentElement).removeClass('text-muted');
            if (failure) {
                $(currentElement).addClass('text-danger');
                $(currentElement).removeClass('empty-circle-muted').addClass('false-circle');

            } else {
                $(currentElement).removeClass('empty-circle-muted').addClass('empty-circle');
            }

            // Steps ahead
        } else {
            $(currentElement).removeClass('text-success').addClass('text-muted');
            $(currentElement).removeClass('check-circle').addClass('empty-circle-muted');

        }

    });
}

// Map step_name to a number for style steps
function display_step(step_name, failure) {
    switch (step_name) {
        case 'submitted':
            style_steps(1, failure);
            break;

        case 'expression_data':
            style_steps(2, failure);
            break;

        case 'ppi_data':
            style_steps(3, failure);
            break;

        case 'validate_preprocess_data':
            style_steps(4, failure);
            break;

        case 'run_clustering':
            style_steps(5, failure);
            break;

        case 'visualize_data':
            style_steps('6', failure);
            break;
    }
}

// Block every alert message except for the one with the id alert_id
function show_only_alert(alert_id) {
    $('.alert').each(function () {
        const alert = $(this);

        // Block every other alert
        if (alert.attr('id') !== alert_id) {
            alert.addClass('d-none')

            // Show selected alert
        } else {
            alert.removeClass('d-none')
        }
    })
}