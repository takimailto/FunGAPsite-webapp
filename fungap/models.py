from django.db import models
import datetime


# Create your models here.


class User(models.Model):
    name = models.CharField(max_length=128, unique=True)
    password = models.CharField(max_length=256)
    email = models.EmailField(unique=True)
    c_time = models.DateTimeField(auto_now_add=True)
    has_confirmed = models.BooleanField(default=False)

    def __str__(self):
        return self.name

    class Meta:
        ordering = ["-c_time"]
        verbose_name = "user"
        verbose_name_plural = "users"


class ConfirmString(models.Model):
    code = models.CharField(max_length=256)
    user = models.OneToOneField('User', on_delete=models.CASCADE)
    c_time = models.DateTimeField(auto_now_add=True)

    def __str__(self):
        return self.user.name + ": " + self.code

    class Meta:
        ordering = ["-c_time"]
        verbose_name = "Confirmation_code"
        verbose_name_plural = "Confirmation_codes"


class Project(models.Model):
    project_name = models.CharField(max_length=128)
    project_description = models.CharField(max_length=1280)
    start_point = models.CharField(max_length=128)
    species = models.CharField(max_length=128)
    u_name = models.ForeignKey(User, on_delete=models.CASCADE)
    status = models.CharField(max_length=128)
    running_information = models.CharField(max_length=1280)  # storage file pathway and so on;

    # json format


def user_directory_path(instance, filename):
    name = instance.project.u_name.name.replace(' ', ' ')
    project_name = instance.project.project_name.replace(' ', '_')
    return 'work_space/users/{0}/{1}/{2}'.format(name, project_name, filename)


class FileUpload(models.Model):
    project = models.ForeignKey(Project, on_delete=models.CASCADE)
    file = models.FileField(upload_to=user_directory_path)


class Custom(models.Model):
    custom_name = models.CharField(max_length=128)
    custom_description = models.CharField(max_length=1280)
    custom_type = models.CharField(max_length=128)
    custom_input_info = models.CharField(max_length=1280)
    custom_output_info = models.CharField(max_length=1280)
    u_name = models.ForeignKey(User, on_delete=models.CASCADE)
