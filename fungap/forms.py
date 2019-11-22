from django import forms
from captcha.fields import CaptchaField
from django.utils import html


class SubmitButtonWidget(forms.Widget):
    def render(self, name, value, attrs=None):
        return '<input type="submit" class="btn btn-primary float-right" name="%s" value="%s">' % (html.escape(name), html.escape(value))


class SubmitButtonField(forms.Field):
    def __init__(self, *args, **kwargs):
        if not kwargs:
            kwargs = {}
        kwargs['widget'] = SubmitButtonWidget

        super(SubmitButtonField, self).__init__(*args, **kwargs)

    def clean(self, value):
        return value


class UserForm(forms.Form):
    username = forms.CharField(label="username", max_length=128, widget=forms.TextInput(attrs={'class': 'form-control',
                                                                                               'placeholder': "Username",
                                                                                               'autofocus': '', 'name': 'uname'}))
    password = forms.CharField(label="password", max_length=256,
                               widget=forms.PasswordInput(attrs={'class': 'form-control', 'placeholder': "Password", 'name': 'psw'}))
    captcha = CaptchaField(label='Verification code')


class RegisterForm(forms.Form):
    username = forms.CharField(label="username", max_length=128, widget=forms.TextInput(attrs={'class': 'form-control',
                                                                                               'placeholder': "Username",
                                                                                               'autofocus': '', 'name': 'uname'}))
    password1 = forms.CharField(label="password", max_length=256,
                                widget=forms.PasswordInput(attrs={'class': 'form-control', 'placeholder': "Password", 'name': 'pws'}))
    password2 = forms.CharField(label="confirm password", max_length=256,
                                widget=forms.PasswordInput(attrs={'class': 'form-control', 'placeholder': "confirm password", 'name': 'pws'}))
    email = forms.EmailField(label="email address", widget=forms.EmailInput(attrs={'class': 'form-comtrol', 'placeholder': 'Email address'}))
    captcha = CaptchaField(label="Verification code")


class ProjectForm(forms.Form):
    projectname = forms.CharField(label="projectname", max_length=128, widget=forms.TextInput(attrs={'class': 'form-control',
                                                                                               'placeholder': "Project name",
                                                                                               'autofocus': ''}))
    projectdescription = forms.CharField(label="projectdescription", max_length=1280, widget=forms.TextInput(attrs={'class': 'form-control',
                                                                                               'placeholder': "Project description",
                                                                                               'autofocus': ''}), required=False)
    CHOICES = (('Fungal', 'Fungal'), ('Other species', 'Other species'),)
    species = forms.ChoiceField(label="species", choices=CHOICES)
    start_CHOICES = (('paired end read', 'paired end read'), ('single end read', 'single end read'), ('bam file', 'bam file'))
    start_point = forms.ChoiceField(label="startpoint", choices=start_CHOICES)


class UploadForm(forms.Form):
    projectname = forms.CharField(label="projectname", max_length=128,
                                  widget=forms.TextInput(attrs={'class': 'form-control',
                                                                'placeholder': "Project name",
                                                                'autofocus': '', "id": "projectnamedefault", "style": "display: none"}))
    paired_end_read = forms.FileField(label="paired_end_read", widget=forms.ClearableFileInput(attrs={'multiple': True}))
    single_end_read = forms.FileField(label="single_end_read")
    bam_file = forms.FileField(label="bam_file")
    augustus_species = forms.CharField(label="augustus_species", max_length=128)
    genome_assembly = forms.FileField(label="genome_assembly")
    sister_proteome = forms.FileField(label="sister_proteome")





