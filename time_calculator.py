import datetime


def time_calculator(start_time, end_time, task):
    total_time = end_time - start_time

    seconds = total_time
    minutes = 0
    hours = 0

    if total_time > 59:
        minutes = total_time // 60
        print(minutes)
        seconds = total_time % 60

        if minutes > 59:
            hours = minutes // 60
            minutes = minutes % 60
    print(f"Time used on {task}: {datetime.time(hours,minutes, seconds)}")
