let error_counter = 0;
const error_threshold = 100;


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

            // Switch depending on status
            switch (task_status) {
                case 'RUNNING':
                    $('#bar-' + task_id).html(progress_percent + '% ');
                    $('#bar-' + task_id).css('width', progress_percent + '%');

                    $('#state-' + task_id).html('RUNNING');
                    setTimeout(update_status, 5000, csrf_token, task_id, poll_url);
                    break;

                case 'SUCCESS':
                    $('#bar-' + task_id).html(100 + '% ');
                    $('#bar-' + task_id).css('width', 100 + '%');
                    $('#bar-' + task_id).addClass('bg-success');

                    $('#state-' + task_id).html('SUCCESS').addClass('text-success');


                    break;

                // Catched exception
                case 'ERROR':
                    $('#bar-' + task_id).html(progress_percent + '% ');
                    $('#bar-' + task_id).css('width', progress_percent + '%');
                    $('#bar-' + task_id).addClass('bg-danger');

                    $('#state-' + task_id).html('ERROR').addClass('text-danger');


                    break;

                // FAILURE means uncatched exception...
                case 'FAILURE':
                    $('#bar-' + task_id).html(0 + '% ');
                    $('#bar-' + task_id).css('width', 0 + '%');
                    $('#bar-' + task_id).addClass('bg-danger');

                    $('#state-' + task_id).html('FAILURE').addClass('text-danger');

                    break;

                case 'PENDING':
                    show_only_alert('alert_pending');
                    $('#state-' + task_id).html('PENDING');
                    setTimeout(update_status, 5000, csrf_token, task_id, poll_url);
                    break;
            }
        },

        error: function (error) {
            if (error_counter > error_threshold){
                return;
            }
            error_counter++;
        }
    });
}