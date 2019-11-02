var poll_xhr;
var finished = 0;
(function () {
    var poll = function () {
        var json_dump = "{{ data }}";
        var task_id = "{{task_id}}";
        console.log(task_id);
        poll_xhr = $.ajax({
            url: 'poll_state',
            type: 'POST',
            data: {
                task_id: task_id,
                csrfmiddlewaretoken: "{{csrf_token}}",
            },
            success: function (result) {
                if (result.process_percent == null || result.process_percent == undefined) {
                    finished = 1;
                    document.getElementById("user-count").textContent = "DONE";
                    jQuery('.bar').css({'width': 100 + '%'});
                    jQuery('.bar').html(100 + '%');
                    document.getElementById('returnBtn').style.visibility = 'visible';
                } else {
                    jQuery('.bar').css({'width': result.process_percent + '%'});
                    jQuery('.bar').html(result.process_percent + '%');
                    document.getElementById("user-count").textContent = "PROCRESSING";
                }
                ;
            }
        });
    };

    var refreshIntervalId = setInterval(function () {
        poll();
        if (finished === 1) {
            clearInterval(refreshIntervalId);
        }
    }, 500);
})();