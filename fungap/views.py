from django.shortcuts import render, redirect, HttpResponse
from fungap import forms, models
import os
import hashlib
import datetime
from django.conf import settings
from django.utils import timezone
from django.contrib import messages
from django.core import serializers
from fungap import tasks
import json
import shutil
import threading
import multiprocessing
from glob import glob
import re

BASE_DIR = '/home/mailto/PycharmProjects/FunGAPsite/work_space/users'


def AnalysisThread(genome_assembly, output_dir, num_cores, trans_read_files, max_intron, sister_proteome, project):
    RunRepeatModeler = multiprocessing.Process(target=tasks.run_repeat_modeler,
                                               args=(genome_assembly, output_dir, num_cores))
    RunHisat2 = multiprocessing.Process(target=tasks.run_hisat2,
                                        args=(genome_assembly, trans_read_files, output_dir, num_cores, max_intron))
    RunHisat2.start()
    RunRepeatModeler.start()
    RunHisat2.join()
    trans_bams2 = []
    for trans_read_file in trans_read_files:
        prefix = re.sub(r'_[12s]$', '',
                        os.path.basename(os.path.splitext(trans_read_file)[0])
                        )
        hisat2_output = os.path.join(os.path.join(output_dir, 'hisat2_out'), '{}.bam'.format(prefix))
        trans_bams2.append(hisat2_output)
    trans_bams = list(set(trans_bams2))
    data = json.loads(project.running_information)
    data['trans_bams'] = 'finished'
    project.running_information = json.dumps(data)
    project.save()
    no_jaccard_clip = '--jaccard_clip'
    RunTrinity = multiprocessing.Process(target=tasks.run_trinity,
                                         args=(trans_bams, output_dir, num_cores, no_jaccard_clip, max_intron))
    RunTrinity.start()
    RunRepeatModeler.join()
    RunTrinity.join()
    repeat_model_file = glob(
        os.path.join(os.path.join(output_dir, 'repeat_modeler_out'), 'RM*/consensi.fa.classified')
    )[0]
    trinity_asms = glob(os.path.join(
        output_dir, 'trinity_out', '*/Trinity_*.fasta')
    )
    data['repeat_model_file'] = 'finished'
    data['trinity_asms'] = 'finished'
    project.running_information = json.dumps(data)
    project.save()
    no_genemark_fungus = '--gmes_fungus'
    augustus_species = json.loads(project.running_information)['augustus_species']
    RunMaker = multiprocessing.Process(target=tasks.run_maker,
                                       args=(genome_assembly, output_dir, augustus_species, sister_proteome, num_cores,
                                             repeat_model_file, trinity_asms, no_genemark_fungus))
    RunMaker.start()
    RunMaker.join()
    masked_assembly = os.path.join(
        output_dir, 'maker_out', 'masked_assembly.fasta'
    )
    RunAugustus = multiprocessing.Process(target=tasks.run_augustus,
                                          args=(masked_assembly, output_dir, augustus_species))
    no_braker_fungus = '--fungus'
    RunBraker1 = multiprocessing.Process(target=tasks.run_braker1,
                                         args=(masked_assembly, trans_bams, output_dir, num_cores, no_braker_fungus))
    RunBraker1.start()
    RunAugustus.start()
    RunBraker1.join()
    RunAugustus.join()
    data['run maker'] = 'finished'
    data['run_augustus'] = 'finished'
    data['run_braker'] = 'finished'
    project.running_information = json.dumps(data)
    project.save()
    maker_faas = glob(os.path.join(output_dir, 'maker_out', '*/maker_*.faa'))
    augustus_faa = os.path.join(output_dir, 'augustus.faa')
    prefixes = [os.path.basename(os.path.splitext(x)[0]) for x in trans_bams]
    prefixes_u = list(set(prefixes))
    braker1_gff3s = []
    braker1_faas = []
    for prefix in prefixes_u:
        braker1_gff3 = os.path.join(
            output_dir, 'braker1_out', prefix, 'braker1_{}.gff3'.format(prefix)
        )
        braker1_gff3s.append(braker1_gff3)
        braker1_faa = os.path.join(
            output_dir, 'braker1_out', prefix, 'braker1_{}.faa'.format(prefix)
        )
        braker1_faas.append(braker1_faa)
    faa_files = [augustus_faa] + maker_faas + braker1_faas
    RunBusco = multiprocessing.Process(target=tasks.run_busco_tatall, args=(faa_files, output_dir, num_cores))
    MakeNrProt = multiprocessing.Process(target=tasks.make_nr_prot, args=(faa_files, output_dir))
    RunBusco.start()
    MakeNrProt.start()
    MakeNrProt.join()
    data['make nr prot'] = 'finished'
    project.running_information = json.dumps(data)
    project.save()
    nr_prot_file = os.path.join(os.path.join(output_dir, 'gene_filtering'), 'nr_prot.faa')
    RunBlastp = multiprocessing.Process(target=tasks.run_blastp,
                                        args=(nr_prot_file, output_dir, sister_proteome, num_cores))
    RunPfamScan = multiprocessing.Process(target=tasks.run_pfam_scan, args=(nr_prot_file, output_dir, num_cores))
    RunBlastp.start()
    RunPfamScan.start()
    trinity_asm = os.path.join(os.path.join(output_dir, 'gene_filtering'), 'trinity_transcripts.fna')
    ConcatenateTransctripts = multiprocessing.Process(target=tasks.concatenate_transcripts,
                                                      args=(trinity_asms, trinity_asm))
    ConcatenateTransctripts.start()
    ConcatenateTransctripts.join()
    maker_gff3s = glob(
        os.path.join(output_dir, 'maker_out', '*/maker_*.gff3')
    )
    augustus_gff3 = os.path.join(output_dir, 'augustus.gff3')
    gff3_files = [augustus_gff3] + maker_gff3s + braker1_gff3s
    RunBlastn = multiprocessing.Process(
        target=tasks.run_blastn_totall(gff3_files, genome_assembly, trinity_asm, output_dir))
    RunBlastn.start()
    BadDit = multiprocessing.Process(target=tasks.catch_bad_genes,
                                     args=(gff3_files, genome_assembly, output_dir))
    BadDit.start()
    RunBlastn.join()
    RunBlastp.join()
    RunPfamScan.join()
    RunBusco.join()
    data['run blastn'] = 'finished'
    data['run pfam scan'] = 'finished'
    data['run blastp'] = 'finished'
    data['run busco'] = 'finished'
    project.running_information = json.dumps(data)
    project.save()
    blastp_output = os.path.join(output_dir, 'gene_filtering', 'nr_prot.blastp')
    busco_out_dir = os.path.join(output_dir, 'busco_out')
    pfam_scan_out = os.path.join(
        output_dir, 'gene_filtering', 'nr_prot.pfam_scan'
    )
    blastn_out_files = []
    for gff3_file in gff3_files:
        blastn_out_files.append('{}.blastn'.format(os.path.join(os.path.join(output_dir, 'gene_filtering'), re.sub(
            r'_transcript\.fna', '', os.path.basename('{}_transcript.fna'.format(os.path.splitext(gff3_file)[0]))
        ))))
    nr_prot_mapping_file = os.path.join(
        os.path.join(output_dir, 'gene_filtering'), 'nr_prot_mapping.txt'
    )
    ImportBlastp = multiprocessing.Process(target=tasks.import_blastp,
                                           args=(blastp_output, nr_prot_mapping_file))
    ImportBusco = multiprocessing.Process(target=tasks.import_busco,
                                          args=(busco_out_dir, output_dir))
    ImportPfam = multiprocessing.Process(target=tasks.import_pfam,
                                         args=(pfam_scan_out, nr_prot_mapping_file))
    ImportBlastn = multiprocessing.Process(target=tasks.import_blastn,
                                           args=(blastn_out_files, output_dir))
    ImportBlastp.start()
    ImportBusco.start()
    ImportPfam.start()
    ImportBlastn.start()
    ImportBlastp.join()
    ImportBusco.join()
    ImportPfam.join()
    ImportBlastn.join()
    BadDit.join()
    data['import blastp'] = 'finished'
    data['import busco'] = 'finished'
    data['import pfam'] = 'finished'
    data['import blastn'] = 'finished'
    data['catch bad gebed'] = 'finished'
    project.running_information = json.dumps(data)
    project.save()
    blastp_dict = os.path.join(os.path.dirname(blastp_output), 'blastp_score.p')
    busco_dict = os.path.join(os.path.join(output_dir, 'gene_filtering'), 'busco_score.p')
    pfam_scan_out_dir = os.path.dirname(pfam_scan_out)
    pfam_dict = os.path.join(pfam_scan_out_dir, 'pfam_score.p')
    blastn_dict = os.path.join(os.path.join(output_dir, 'gene_filtering'), 'blastn_score.p')
    bad_dict = os.path.join(os.path.join(output_dir, 'gene_filtering'), 'D_bad.p')
    tasks.filter_gff3s(
        genome_assembly, gff3_files, blastp_dict, busco_dict, pfam_dict,
        blastn_dict, bad_dict, nr_prot_file, nr_prot_mapping_file, output_dir
    )
    tasks.gff3_postprocess(genome_assembly, os.path.join(output_dir, 'gene_filtering', 'filtered_1.gff3'))
    # Copy output files
    tasks.copy_output(output_dir)
    # Create markdown
    tasks.create_markdown(genome_assembly, output_dir, trans_bams, trinity_asms, '')
    project.status = 'finished'
    project.save()
    os.system('mkdir {}'.format(BASE_DIR))
    static_path = os.path.join('/home/mailto/PycharmProjects/FunGAPsite/static/fungap/images', project.project_name)
    if not os.path.exists(static_path):
        os.system('mkdir {}'.format(static_path))
    ori_path = os.path.join(BASE_DIR, 'Ace', project.project_name, 'output', 'fungap_out')
    os.system('cp {} {}'.format(os.path.join(ori_path, 'fungap_out_prot_len_dist.png'), static_path))
    os.system('cp {} {}'.format(os.path.join(ori_path, 'fingap_out_trans_len_disr.png'), static_path))


def homepage(request):
    login=True
    upload_form = forms.UploadForm()
    login_form = forms.UserForm()
    register_form = forms.RegisterForm()
    if request.method == "POST" and "signin" in request.POST:
        if request.session.get('is_login', None):  # Do not allow duplicate login\
            message = "You have sign in, do not duplicate login"
            return render(request, 'fungap/homepage.html', locals())
        else:
            login_form = forms.UserForm(request.POST)
            message = 'Please check the completed content!'
            if login_form.is_valid():
                username = login_form.cleaned_data.get('username')
                password = login_form.cleaned_data.get('password')
                try:
                    user = models.User.objects.get(name=username)
                except:
                    message = 'Username does not exist!'
                    return render(request, 'fungap/homepage.html', locals())

                if not user.has_confirmed:
                    message = 'This user has not been confirmed by email!'
                    return render(request, 'fungap/homepage.html', locals())

                if user.password == hash_code(password):
                    request.session['is_login'] = True
                    request.session['user_id'] = user.id
                    request.session['user_name'] = user.name

                    return redirect('/homepage/')
                else:
                    message = 'The password is incorrect!'
                    return render(request, 'fungap/homepage.html', locals())  # locals() :Return all current local
                    # variable dictionaries
            else:
                return render(request, 'fungap/homepage.html', locals())
    if request.method == "POST" and "register" in request.POST:
        if request.session.get('is_login', None):
            return render(request, 'fungap/homepage.html', locals())
        else:
            register_form = forms.RegisterForm(request.POST)
            message = "Please check the completed content!"
            if register_form.is_valid():
                username = register_form.cleaned_data.get('username')
                password1 = register_form.cleaned_data.get('password1')
                password2 = register_form.cleaned_data.get('password2')
                email = register_form.cleaned_data.get('email')

                if password1 != password2:
                    message = 'The password entered twice is different!'
                    return render(request, 'fungap/homepage.html', locals())
                else:
                    same_name_usr = models.User.objects.filter(name=username)
                    if same_name_usr:
                        message = 'Username already exists!'
                        return render(request, 'fungap/homepage.html', locals())
                    same_email_user = models.User.objects.filter(email=email)
                    if same_email_user:
                        message = 'The mailbox has been registered!'
                        return render(request, 'fungap/homepage.html', locals())

                    new_user = models.User()
                    new_user.name = username
                    new_user.password = hash_code(password1)
                    new_user.email = email
                    new_user.save()

                    code = make_confirm_string(new_user)
                    send_email(email, code)

                    message = 'Please go to the email to confirm!'
                    return render(request, 'fungap/confirm.html', locals())
            else:
                return render(request, 'fungap/homepage.html', locals())
    if request.method == 'POST' and 'signout' in request.POST:
        if not request.session.get('is_login', None):
            return render(request, 'fungap/homepage.html', locals())
        request.session.flush()
        # or
        # del request.session['is_login']
        # del request.session['user_id']
        # del request.session['user_name']
        return render(request, "fungap/homepage.html", locals())
    if request.method == "POST" and 'upload' in request.POST:
        uploadform = forms.UploadForm(request.POST, request.FILES)
        uploadform.is_valid()
        projects = models.User.objects.get(name='Ace').project_set.all()
        project = projects.get(project_name=uploadform.cleaned_data.get('projectname'))
        DIR = os.path.join(BASE_DIR, project.u_name.name.replace(' ', '_'), project.project_name.replace(' ', '_'))
        if project.start_point == 'paired end read':
            trans_read_single = ''
            trans_bam = ''
            files = request.FILES.getlist('paired_end_read')
            paired_end_read = []
            for file in files:
                fileupload = models.FileUpload(project=project, file=file)
                fileupload.save()
                fileupload.delete()
                path = os.path.join(DIR, file.name)
                paired_end_read.append(os.path.join(DIR, file.name))
            if len(paired_end_read) != 2:
                message = "Please upload 2 paired end reads!"
                shutil.rmtree(DIR)
                messages.add_message(request, messages.INFO,
                                     message)
                return redirect('/homepage/')
            trans_read_1 = paired_end_read[0]
            trans_read_2 = paired_end_read[1]
        elif project.start_point == 'single end read':
            fileupload = models.FileUpload(project=project, file=request.FILES['single_end_read'])
            fileupload.save()
            fileupload.delete()
            trans_read_1 = ''
            trans_read_2 = ''
            trans_bam = ''
            trans_read_single = os.path.join(DIR, request.FILES['single_end_read'].name)
        else:
            fileupload = models.FileUpload(project=project, file=request.FILES['bam_file'])
            fileupload.save()
            fileupload.delete()
            trans_read_1 = ''
            trans_read_2 = ''
            trans_read_single = ''
            trans_bam = os.path.join(DIR, request.FILES['bam_file'].name)
        fileupload = models.FileUpload(project=project, file=request.FILES['genome_assembly'])
        fileupload.save()
        fileupload.delete()
        genome_assembly = os.path.join(DIR, request.FILES['genome_assembly'].name)
        fileupload = models.FileUpload(project=project, file=request.FILES['sister_proteome'])
        fileupload.save()
        fileupload.delete()
        sister_proteome = os.path.join(DIR, request.FILES['sister_proteome'].name)
        ori = json.loads('{}')
        ori['augustus_species'] = request.POST['augustus_species']
        project.running_information = json.dumps(ori)
        project.save()
        tasks.set_loggings(os.path.join(DIR, 'output'))
        trans_read_files = tasks.check_inputs(trans_read_1, trans_read_2, trans_read_single, trans_bam, genome_assembly,
                                              sister_proteome)
        if type(trans_read_files) == str:
            message = trans_read_files + ' Please upload again.'
            shutil.rmtree(DIR)
            messages.add_message(request, messages.INFO,
                                 message)
            return redirect('/homepage/')
        else:
            message = 'FunGAP is running now. Wait...'
            messages.add_message(request, messages.INFO, message)
            project.status = 'Running'
            project.save()
            analysis = threading.Thread(target=AnalysisThread,
                                        args=(genome_assembly, os.path.join(DIR, 'output'), 1, trans_read_files, 2000,
                                              sister_proteome, project))
            analysis.start()
            return redirect('/homepage/')
    else:
        data = serializers.serialize("json",
                                     models.User.objects.get(name='Ace').project_set.all())
        return render(request, 'fungap/homepage.html', locals())


# Password encryption


def hash_code(s, salt='FunGAPsite'):
    h = hashlib.sha256()
    s += salt
    h.update(s.encode())
    return h.hexdigest()


# Create a confirmation code object


def make_confirm_string(user):
    now = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    code = hash_code(user.name, now)
    models.ConfirmString.objects.create(code=code, user=user)
    return code


def send_email(email, code):
    from django.core.mail import EmailMultiAlternatives

    subject = "Registration confirmation email from FunGAPsite."
    text_content = '''Thanks for register FunGAPsite. fGAP (fungal Genome Annotation Pipeline) predicts gene features 
    from genome assembly and RNA-seq reads. Current version only aims to fungal genome, but we will extend the scope 
    to larger eukaryotic genomes later. '''
    html_content = '''<p>Thanks for register <a href="http://{}/confirm/?code={}" target=blank>FunGAPsite</a>, 
    fGAP (fungal Genome Annotation Pipeline) predicts gene features from genome assembly and RNA-seq reads. Current 
    version only aims to fungal genome, but we will extend the scope to larger eukaryotic genomes later.</p> 
    <p>Please click on the site link to complete the registration confirmation!</p> <p>This link is valid for {} 
    days!</p> '''.format('127.0.0.1:8000', code, settings.CONFIRM_DAYS)

    msg = EmailMultiAlternatives(subject, text_content, settings.EMAIL_HOST_USER, [email])
    msg.attach_alternative(html_content, "text/html")
    msg.send()


def user_confirm(request):
    code = request.GET.get('code', None)
    message = ''
    try:
        confirm = models.ConfirmString.objects.get(code=code)
    except:
        message = 'Invalid confirmation request!'
        return render(request, 'fungap/confirm.html', locals())

    c_time = confirm.c_time
    now = timezone.now()
    if now > c_time + datetime.timedelta(settings.CONFIRM_DAYS):
        confirm.user.delete()
        message = 'Your email has expired! Please re-register!'
        return render(request, 'fungap/confirm.html', locals())
    else:
        confirm.user.has_confirmed = True
        confirm.user.save()
        confirm.delete()
        message = 'Thanks for confirming, please log in with your account!'
        return render(request, 'fungap/confirm.html', locals())


def new(request):
    project_form = forms.ProjectForm
    if request.method == "POST" and 'create' in request.POST:
        project_form = forms.ProjectForm(request.POST)
        message = 'Please check the completed content!'
        if project_form.is_valid():
            projectname = project_form.cleaned_data.get('projectname')
            projectdescription = project_form.cleaned_data.get('projectdescription')
            species = project_form.cleaned_data.get('species')
            startpoint = project_form.cleaned_data.get('start_point')
            projects = models.User.objects.get(name='Ace').project_set.all()
            for project in projects:
                if projectname == project.project_name:
                    message = "Project name already exists"
                    return render(request, 'fungap/new.html', locals())
            new_project = models.Project()
            new_project.project_name = projectname
            new_project.project_description = projectdescription
            new_project.start_point = startpoint
            new_project.species = species
            new_project.u_name = models.User.objects.get(name='Ace')
            new_project.status = 'upload file'
            new_project.save()
            messages.add_message(request, messages.INFO,
                                 'Project has been created successfully, please upload your work now.')
            return redirect('/homepage/')
        else:
            return render(request, 'fungap/', locals())  # locals() :Return all current local
    if request.method == 'POST' and 'signout' in request.POST:
        if not request.session.get('is_login', None):
            return redirect('/homepage/')
        request.session.flush()
        # or
        # del request.session['is_login']
        # del request.session['user_id']
        # del request.session['user_name']
        return redirect('/homepage/')
    return render(request, 'fungap/new.html', locals())


def custom(request):
    if request.method == "POST" and 'hisat2' in request.POST:
        custom_reads = request.POST['reads'].split(' ')
        custom_max_intron = request.POST['maxintron']
        custom_genome_assembly = request.POST['genomeassembly']
        custom_output_dir = request.POST['outputdir']
        new_custom = models.Custom()
        new_custom.custom_type = 'hisat2'
        data = {}
        data1 = {}
        data['reads'] = request.POST['reads']
        data['max_intron'] = custom_max_intron
        data['genome_assembly'] = custom_genome_assembly
        data1['output_dir'] = custom_output_dir
        data1['running_status'] = 'start'
        new_custom.custom_output_info = json.dumps(data1)
        new_custom.custom_input_info = json.dumps(data)
        new_custom.u_name = models.User.objects.get(name='Ace')
        new_custom.save()
        messages.add_message(request, messages.INFO,
                             'Custom hista2 project is running now.')
        hisat2 = threading.Thread(target=hisat2Thread, args=(
            custom_genome_assembly, custom_reads, custom_output_dir, 1, custom_max_intron,
            new_custom.id))
        hisat2.start()
        return redirect('/custom/')
    if request.method == 'POST' and 'trinity_run' in request.POST:
        custom_bams = request.POST['bams'].split(' ')
        custom_max_intron = request.POST['maxintron']
        custom_output_dir = request.POST['outputdir']
        new_custom = models.Custom()
        new_custom.custom_type = 'trinity'
        data = {}
        data1 = {}
        data['bams'] = request.POST['bams']
        data['max_intron'] = custom_max_intron
        data1['output_dir'] = custom_output_dir
        data1['running_status'] = 'start'
        if request.POST.get('jaccardclip'):
            no_jaccard_clip = '--jaccard_clip'
            data['jaccard clip'] = 'yes'
        else:
            no_jaccard_clip = ''
            data['jaccard clip'] = 'no'
        new_custom.custom_output_info = json.dumps(data1)
        new_custom.custom_input_info = json.dumps(data)
        new_custom.u_name = models.User.objects.get(name="Ace")
        new_custom.save()
        messages.add_message(request, messages.INFO,
                             'Custom Trinity project is running now.')
        trinity = threading.Thread(target=trinityThread, args=(
            custom_bams, custom_output_dir, 1, no_jaccard_clip, custom_max_intron, new_custom.id))
        trinity.start()
        return redirect('/custom/')
    if request.method == 'POST' and 'repeatmodeler' in request.POST:
        custom_genome_assembly = request.POST['genomeassembly']
        custom_output_dir = request.POST['outputdir']
        new_custom = models.Custom()
        new_custom.custom_type = 'repeat modeler'
        data = {}
        data1 = {}
        data['genome_assembly'] = custom_genome_assembly
        data1['output_dir'] = custom_output_dir
        data1['running_status'] = 'start'
        new_custom.custom_output_info = json.dumps(data1)
        new_custom.custom_input_info = json.dumps(data)
        new_custom.u_name = models.User.objects.get(name='Ace')
        new_custom.save()
        messages.add_message(request, messages.INFO,
                             'Custom repeat modeler project is running now.')
        repeat_modeler = threading.Thread(target=RepeatModelerThread, args=(
            custom_genome_assembly, custom_output_dir, 1, new_custom.id))
        repeat_modeler.start()
        return redirect('/custom/')
    if request.method == "POST" and 'repeatmasker' in request.POST:
        custom_genome_assembly = request.POST['genomeassembly']
        custom_repeat_model_file = request.POST['repeatmodel']
        custom_output_dir = request.POST['outputdir']
        new_custom = models.Custom()
        new_custom.custom_type = 'RepeatMasker'
        data = {}
        data1 = {}
        data['genome_assembly'] = custom_genome_assembly
        data['repeat_model_file'] = custom_repeat_model_file
        data1['output_dir'] = custom_output_dir
        data1['running_status'] = 'start'
        new_custom.custom_output_info = json.dumps(data1)
        new_custom.custom_input_info = json.dumps(data)
        new_custom.u_name = models.User.objects.get(name='Ace')
        new_custom.save()
        messages.add_message(request, messages.INFO,
                             'Custom repeat masker project is running now.')
        RepeatMasker = threading.Thread(target=repeatmaskerThread, args=(
            custom_genome_assembly, custom_output_dir, custom_repeat_model_file, 1, new_custom.id))
        RepeatMasker.start()
        return redirect('/custom/')
    if request.method == 'POST' and 'maker' in request.POST:
        custom_genome_assembly = request.POST['genomeassembly']
        custom_repeat_model_file = request.POST['repeatmodel']
        custom_augustus_species = request.POST['augustus']
        custom_sister_proteome = request.POST['sisterprotein']
        custom_trinity_asms = request.POST['trinityasms'].split(' ')
        custom_output_dir = request.POST['outputdir']
        new_custom = models.Custom()
        new_custom.custom_type = 'maker'
        data = {}
        data1 = {}
        data['genome_assembly'] = custom_genome_assembly
        data['repeat_model_file'] = custom_repeat_model_file
        data['augustus_species'] = custom_augustus_species
        data['sister_proteom'] = custom_sister_proteome
        data['trinity_asms'] = request.POST['trinityasms']
        data1['output_dir'] = custom_output_dir
        data1['running_status'] = 'start'
        if request.POST.get('genemarkfungus'):
            no_genemark_fungus = '--gmes_fungus'
            data['genemark fungus'] = 'yes'
        else:
            no_genemark_fungus = ''
            data['genemark fungus'] = 'no'
        new_custom.custom_output_info = json.dumps(data1)
        new_custom.custom_input_info = json.dumps(data)
        new_custom.u_name = models.User.objects.get(name="Ace")
        new_custom.save()
        messages.add_message(request, messages.INFO,
                             'Custom maker project is running now.')
        maker = threading.Thread(target=makerThread, args=(
            custom_genome_assembly, custom_output_dir, custom_augustus_species, custom_sister_proteome, 4,
            custom_repeat_model_file, custom_trinity_asms, no_genemark_fungus, new_custom.id))
        maker.start()
        return redirect('/custom/')
    if request.method == "POST" and 'braker1' in request.POST:
        custom_masked_assembly = request.POST['maskedassembly']
        custom_bams = request.POST['bams'].split(' ')
        custom_output_dir = request.POST['outputdir']
        new_custom = models.Custom()
        new_custom.custom_type = 'braker1'
        data = {}
        data1 = {}
        data['masked_assembly'] = custom_masked_assembly
        data['bams'] = request.POST['bams']
        data1['output_dir'] = custom_output_dir
        data1['running_status'] = 'start'
        if request.POST.get('brakerfungus'):
            no_braker_fungus = '--fungus'
            data['braker fungus'] = 'yes'
        else:
            no_braker_fungus = ''
            data['braker fungus'] = 'no'
        new_custom.custom_output_info = json.dumps(data1)
        new_custom.custom_input_info = json.dumps(data)
        new_custom.u_name = models.User.objects.get(name="Ace")
        new_custom.save()
        messages.add_message(request, messages.INFO,
                             'Custom braker1 project is running now.')
        maker = threading.Thread(target=brakerThread, args=(
            custom_masked_assembly, custom_bams, custom_output_dir, 1, no_braker_fungus, new_custom.id))
        maker.start()
        return redirect('/custom/')
    if request.method == "POST" and 'augustus_run' in request.POST:
        custom_masked_assembly = request.POST['maskedassembly']
        custom_augustus_species = request.POST['augustus']
        custom_output_dir = request.POST['outputdir']
        new_custom = models.Custom()
        new_custom.custom_type = 'augustus'
        data = {}
        data1 = {}
        data['masked_assembly'] = custom_masked_assembly
        data['augustus_species'] = custom_augustus_species
        data1['output_dir'] = custom_output_dir
        data1['running_status'] = 'start'
        new_custom.custom_output_info = json.dumps(data1)
        new_custom.custom_input_info = json.dumps(data)
        new_custom.u_name = models.User.objects.get(name="Ace")
        new_custom.save()
        messages.add_message(request, messages.INFO,
                             'Custom augustus project is running now.')
        augustus = threading.Thread(target=augustusThread, args=(
            custom_masked_assembly, custom_output_dir, custom_augustus_species, new_custom.id))
        augustus.start()
        return redirect('/custom/')
    if request.method == "POST" and 'busco' in request.POST:
        custom_faas = request.POST['faas'].split(' ')
        custom_output_dir = request.POST['outputdir']
        new_custom = models.Custom()
        new_custom.custom_type = 'busco'
        data = {}
        data1 = {}
        data['faas'] = request.POST['faas']
        data1['output_dir'] = custom_output_dir
        data1['running_status'] = 'start'
        new_custom.custom_output_info = json.dumps(data1)
        new_custom.custom_input_info = json.dumps(data)
        new_custom.u_name = models.User.objects.get(name="Ace")
        new_custom.save()
        messages.add_message(request, messages.INFO,
                             'Custom busco project is running now.')
        busco = threading.Thread(target=buscoThread, args=(
            custom_faas, custom_output_dir, 1, new_custom.id))
        busco.start()
        return redirect('/custom/')
    if request.method == "POST" and 'nr_run' in request.POST:
        custom_faas = request.POST['faas'].split(' ')
        custom_output_dir = request.POST['outputdir']
        new_custom = models.Custom()
        new_custom.custom_type = 'make nr proteins'
        data = {}
        data1 = {}
        data['faas'] = request.POST['faas']
        data1['output_dir'] = custom_output_dir
        data1['running_status'] = 'start'
        new_custom.custom_output_info = json.dumps(data1)
        new_custom.custom_input_info = json.dumps(data)
        new_custom.u_name = models.User.objects.get(name="Ace")
        new_custom.save()
        messages.add_message(request, messages.INFO,
                             'Custom Create none redundant protein database project is running now.')
        nr = threading.Thread(target=nrThread, args=(
            custom_faas, custom_output_dir, new_custom.id))
        nr.start()
        return redirect('/custom/')
    if request.method == "POST" and 'blastp' in request.POST:
        custom_nr_prot = request.POST['nrprot']
        custom_output_dir = request.POST['outputdir']
        custom_sister_proteome = request.POST['sisterprotein']
        new_custom = models.Custom()
        new_custom.custom_type = 'blastp'
        data = {}
        data1 = {}
        data['sister_proteom'] = custom_sister_proteome
        data['nr_prot'] = custom_nr_prot
        data1['output_dir'] = custom_output_dir
        data1['running_status'] = 'start'
        new_custom.custom_output_info = json.dumps(data1)
        new_custom.custom_input_info = json.dumps(data)
        new_custom.u_name = models.User.objects.get(name="Ace")
        new_custom.save()
        messages.add_message(request, messages.INFO,
                             'Custom Blastp project is running now.')
        blastp = threading.Thread(target=blastpThread, args=(
            custom_nr_prot, custom_output_dir, custom_sister_proteome, new_custom.id))
        blastp.start()
        return redirect('/custom/')
    if request.method == "POST" and 'interpro' in request.POST:
        custom_nr_prot = request.POST['nrprot']
        custom_output_dir = request.POST['outputdir']
        new_custom = models.Custom()
        new_custom.custom_type = 'InterProScan'
        data = {}
        data1 = {}
        data['nr_prot'] = custom_nr_prot
        data1['output_dir'] = custom_output_dir
        data1['running_status'] = 'start'
        new_custom.custom_output_info = json.dumps(data1)
        new_custom.custom_input_info = json.dumps(data)
        new_custom.u_name = models.User.objects.get(name="Ace")
        new_custom.save()
        messages.add_message(request, messages.INFO,
                             'Custom pfam scan project is running now.')
        pfamscan = threading.Thread(target=pfamscanThread, args=(
            custom_nr_prot, custom_output_dir, new_custom.id))
        pfamscan.start()
        return redirect('/custom/')
    if request.method == 'POST' and 'blastn' in request.POST:
        custom_genome_assembly = request.POST['genomeassembly']
        custom_gff3s = request.POST['gff3s'].split(' ')
        custom_trinity_asm = request.POST['trinity_asm']
        custom_output_dir = request.POST['outputdir']
        new_custom = models.Custom()
        new_custom.custom_type = 'blastn'
        data = {}
        data1 = {}
        data['genome_assembly'] = custom_genome_assembly
        data1['output_dir'] = custom_output_dir
        data['gff3s'] = request.POST['gff3s']
        data['trinity_asm'] = custom_trinity_asm
        data1['running_status'] = 'start'
        new_custom.custom_output_info = json.dumps(data1)
        new_custom.custom_input_info = json.dumps(data)
        new_custom.u_name = models.User.objects.get(name="Ace")
        new_custom.save()
        messages.add_message(request, messages.INFO,
                             'Custom blatsn project is running now.')
        blastn = threading.Thread(target=blastnThread, args=(
            custom_gff3s, custom_genome_assembly, custom_trinity_asm, custom_output_dir, new_custom.id))
        blastn.start()
        return redirect('/custom/')
    if request.method == 'POST' and 'concatenate' in request.POST:
        trinity_asms = request.POST['trinityasms'].split(' ')
        trinity_asm = os.path.join(os.path.dirname(trinity_asms[0]), 'gene_filtering', 'trinity_transcripts.fna')
        concatenate = threading.Thread(target=tasks.concatenate_transcripts, args=(trinity_asms, trinity_asm))
        concatenate.start()
        messages.add_message(request, messages.INFO,
                             'Finished, the result is in {}.'.format(trinity_asm))
        return redirect('/custom/')
    if request.method == "POST" and 'importbp' in request.POST:
        custom_nr_mapping = request.POST['nrprotmapping']
        custom_blastp_out = request.POST['importbp_bp']
        new_custom = models.Custom()
        data = {}
        data1 = {}
        data['nr_prot_mapping'] = custom_nr_mapping
        data['blastp'] = custom_blastp_out
        data1['running_status'] = 'start'
        new_custom.custom_output_info = json.dumps(data1)
        new_custom.custom_input_info = json.dumps(data)
        new_custom.u_name = models.User.objects.get(name="Ace")
        new_custom.custom_type = 'import blastp'
        new_custom.custom_name = ' '
        new_custom.save()
        messages.add_message(request, messages.INFO,
                             'Custom import_blastp project is running now.')
        import_blastp = threading.Thread(target=importblastpThread,
                                         args=(custom_blastp_out, custom_nr_mapping, new_custom.id))
        import_blastp.start()
        return redirect('/custom/')
    if request.method == "POST" and 'importbo' in request.POST:
        custom_busco_out = request.POST['importbo_bo']
        new_custom = models.Custom()
        data = {}
        data1 = {}
        data['busco'] = custom_busco_out
        data1['running_status'] = 'start'
        new_custom.custom_output_info = json.dumps(data1)
        new_custom.custom_input_info = json.dumps(data)
        new_custom.u_name = models.User.objects.get(name="Ace")
        new_custom.custom_type = 'import busco'
        new_custom.custom_name = ' '
        new_custom.save()
        messages.add_message(request, messages.INFO,
                             'Custom import_blastp project is running now.')
        import_busco = threading.Thread(target=importbuscoThread, args=(custom_busco_out, new_custom.id))
        import_busco.start()
        return redirect('/custom/')
    if request.method == "POST" and 'importpm' in request.POST:
        custom_nr_mapping = request.POST['nrprotmapping']
        custom_pfam_out = request.POST['importpm_pm']
        new_custom = models.Custom()
        data = {}
        data1 = {}
        data['nr_prot_mapping'] = custom_nr_mapping
        data['pfam'] = custom_pfam_out
        data1['running_status'] = 'start'
        new_custom.custom_output_info = json.dumps(data1)
        new_custom.custom_input_info = json.dumps(data)
        new_custom.u_name = models.User.objects.get(name="Ace")
        new_custom.custom_type = 'import pfam scan'
        new_custom.custom_name = ' '
        new_custom.save()
        messages.add_message(request, messages.INFO,
                             'Custom import_blastp project is running now.')
        import_pfam = threading.Thread(target=importpfamThread,
                                       args=(custom_pfam_out, custom_nr_mapping, new_custom.id))
        import_pfam.start()
        return redirect('/custom/')
    if request.method == "POST" and 'importbn' in request.POST:
        custom_blastn_output_files = request.POST['importbn_bn'].split(' ')
        new_custom = models.Custom()
        data = {}
        data1 = {}
        data['custom_blastn_output_file'] = request.POST['importbn_bn']
        data1['running_status'] = 'start'
        new_custom.custom_output_info = json.dumps(data1)
        new_custom.custom_input_info = json.dumps(data)
        new_custom.u_name = models.User.objects.get(name="Ace")
        new_custom.custom_type = 'import pfam scan'
        new_custom.custom_name = ' '
        new_custom.save()
        messages.add_message(request, messages.INFO,
                             'Custom import_blastp project is running now.')
        import_blastn = threading.Thread(target=importblastnThread,
                                         args=(custom_blastn_output_files, new_custom.id))
        import_blastn.start()
        return redirect('/custom/')
    if request.method == "POST" and 'catchbad' in request.POST:
        custom_output_dir = request.POST['outputdir']
        custom_genome_assembly = request.POST['genomeassembly']
        custom_gff3s = request.POST['gff3s'].split(' ')
        data = {}
        data1 = {}
        data['genome_assembly'] = custom_genome_assembly
        data['gff3s'] = request.POST['gff3s']
        data1['running_status'] = 'start'
        data1['output_dir'] = custom_output_dir
        new_custom = models.Custom()
        new_custom.custom_output_info = json.dumps(data1)
        new_custom.custom_input_info = json.dumps(data)
        new_custom.u_name = models.User.objects.get(name="Ace")
        new_custom.custom_type = 'cathch bad genes'
        new_custom.custom_name = ' '
        new_custom.save()
        messages.add_message(request, messages.INFO,
                             'Custom catch bad genes project is running now.')
        catch_bad_genes = threading.Thread(target=baddictThread,
                                           args=(
                                               custom_gff3s, custom_genome_assembly, custom_output_dir, new_custom.id))
        catch_bad_genes.start()
        return redirect('/custom/')
    if request.method == "POST" and 'filter' in request.POST:
        custom_nr_mapping = request.POST['noneredundantprotmapping']
        custom_nr_prot = request.POST['noneredundantprot']
        custom_bad_dict = request.POST['bdict']
        custom_blastn_dict = request.POST['bndict']
        custom_blastp_dict = request.POST['bpdict']
        custom_pfam_dict = request.POST['pmdict']
        custom_busco_dict = request.POST['bodict']
        custom_output_dir = request.POST['outputdir']
        custom_genome_assembly = request.POST['genomeassembly']
        custom_gff3s = request.POST['gff3s'].split(' ')
        new_custom = models.Custom()
        data = {}
        data1 = {}
        data['nr_prot_mapping'] = custom_nr_mapping
        data['nr_prot'] = custom_nr_prot
        data['bad_dict'] = custom_bad_dict
        data['blastn_dict'] = custom_blastn_dict
        data['blastp_dict'] = custom_blastp_dict
        data['pfam_dict'] = custom_pfam_dict
        data['busco_dict'] = custom_busco_dict
        data['genome_assembly'] = custom_genome_assembly
        data['gff3s'] = request.POST['gff3s']
        data1['running_status'] = 'start'
        data1['output_dir'] = custom_output_dir
        new_custom.custom_output_info = json.dumps(data1)
        new_custom.custom_input_info = json.dumps(data)
        new_custom.u_name = models.User.objects.get(name="Ace")
        new_custom.custom_type = 'filter'
        new_custom.custom_name = ' '
        new_custom.save()
        messages.add_message(request, messages.INFO,
                             'Custom filter project is running now.')
        filter = threading.Thread(target=filterThread, args=(
            custom_genome_assembly, custom_gff3s, custom_blastp_dict, custom_busco_dict, custom_pfam_dict,
            custom_blastn_dict, custom_bad_dict, custom_nr_prot, custom_nr_mapping, custom_output_dir, new_custom.id))
        filter.start()
        return redirect('/custom/')
    if request.method == "POST" and 'postprocess' in request.POST:
        custom_genome_assemply = request.POST['genomeassembly']
        custom_ori_gff3 = request.POST['ori_gff3']
        new_custom = models.Custom()
        data = {}
        data1 = {}
        data['ori_gff'] = custom_ori_gff3
        data['genome_assembly'] = custom_genome_assemply
        data1['running_status'] = 'start'
        new_custom.custom_output_info = json.dumps(data1)
        new_custom.custom_input_info = json.dumps(data)
        new_custom.u_name = models.User.objects.get(name="Ace")
        new_custom.custom_type = 'post process'
        new_custom.custom_name = ' '
        new_custom.save()
        messages.add_message(request, messages.INFO,
                             'Custom post process project is running now.')
        post_process = threading.Thread(target=postprocessThread,
                                        args=(custom_genome_assemply, custom_ori_gff3, new_custom.id))
        post_process.start()
        return redirect('/custom/')
    if request.method == "POST" and 'createmarkdown' in request.POST:
        custom_genome_assemply = request.POST['genomeassembly']
        custom_final_gff3 = request.POST['final_gff3']
        custom_trinity_asms = request.POST['trinity_asms'].split(' ')
        custom_bam_files = request.POST['bams'].split(' ')
        custom_output_dir = request.POST['outputdir']
        new_custom = models.Custom()
        data = {}
        data1 = {}
        data['final_gff'] = custom_final_gff3
        data['genome_assembly'] = custom_genome_assemply
        data['trinity_asms'] = request.POST['trinity_asms']
        data['bam_files'] = request.POST['bams']
        data1['running_status'] = 'start'
        new_custom.custom_output_info = json.dumps(data1)
        new_custom.custom_input_info = json.dumps(data)
        new_custom.u_name = models.User.objects.get(name="Ace")
        new_custom.custom_type = 'create markdown'
        new_custom.custom_name = ' '
        new_custom.save()
        messages.add_message(request, messages.INFO,
                             'Custom create markdown process project is running now.')
        create_markdown = threading.Thread(target=createmarkdownThread, args=(
            custom_genome_assemply, custom_final_gff3, custom_trinity_asms, custom_bam_files, custom_output_dir,
            new_custom.id))
        create_markdown.start()
        return redirect('/custom/')
    if request.method == "POST" and "downloadprotdb" in request.POST:
        taxon = request.POST['taxon']
        email_address = request.POST['email']
        output_dir = request.POST['outputdir']
        num_sisters = request.POST['num_sisters']
        download_prot_db = threading.Thread(target=tasks.download_sister_orgs,
                                            args=(taxon, email_address, num_sisters, output_dir))
        download_prot_db.start()
        messages.add_message(request, messages.INFO,
                             'Download successfully, Go to your output dir to check it!')
        return redirect('/custom/')
    if request.method == "POST" and "getaugustuspecies" in request.POST:
        genus_name = request.POST['genusname']
        email_address = request.POST['email']
        get_augustus_species = threading.Thread(target=tasks.get_augustus_specie, args=(genus_name, email_address))
        get_augustus_species.start()
        messages.add_message(request, messages.INFO,
                             'You can see the augustus species in the command line.')
        return redirect('/custom/')
    if request.method == 'POST' and 'signout' in request.POST:
        if not request.session.get('is_login', None):
            return redirect('/homepage/')
        request.session.flush()
        # or
        # del request.session['is_login']
        # del request.session['user_id']
        # del request.session['user_name']
        return redirect('/homepage/')
    else:
        running_info = serializers.serialize("json",
                                             models.User.objects.get(
                                                 name="Ace").custom_set.all())
        return render(request, 'fungap/custom.html', locals())


def hisat2Thread(genome_assembly, trans_read_files, output_dir, num_cores, max_intron, custom_id):
    custom = models.Custom.objects.get(id=custom_id)
    tasks.set_loggings(os.path.join(output_dir, 'output'))
    tasks.run_hisat2(genome_assembly, trans_read_files, output_dir, num_cores, max_intron)
    trans_bams2 = []
    for trans_read_file in trans_read_files:
        prefix = re.sub(r'_[12s]$', '',
                        os.path.basename(os.path.splitext(trans_read_file)[0])
                        )
        hisat2_output = os.path.join(os.path.join(output_dir, 'hisat2_out'), '{}.bam'.format(prefix))
        trans_bams2.append(hisat2_output)
    trans_bams = list(set(trans_bams2))
    data1 = {}
    data1['running_status'] = 'finished'
    data1['output_files'] = trans_bams
    custom.custom_output_info = json.dumps(data1)
    custom.save()


def trinityThread(trans_bams, output_dir, num_cores, no_jaccard_clip, max_intron, custom_id):
    custom = models.Custom.objects.get(id=custom_id)
    tasks.set_loggings(os.path.join(output_dir, 'output'))
    tasks.run_trinity(trans_bams, output_dir, num_cores, no_jaccard_clip, max_intron)
    trinity_asms = glob(os.path.join(
        output_dir, 'trinity_out', '*/Trinity_*.fasta')
    )
    data1 = {}
    data1['running_status'] = 'finished'
    data1['output_files'] = trinity_asms
    custom.custom_output_info = json.dumps(data1)
    custom.save()


def RepeatModelerThread(genome_assembly, output_dir, num_cores, custom_id):
    custom = models.Custom.objects.get(id=custom_id)
    tasks.set_loggings(os.path.join(output_dir, 'output'))
    tasks.run_repeat_modeler(genome_assembly, output_dir, num_cores)
    repeat_model_file = glob(
        os.path.join(os.path.join(output_dir, 'repeat_modeler_out'), 'RM*/consensi.fa.classified')
    )[0]
    data1 = {}
    data1['running_status'] = 'finished'
    data1['output_files'] = repeat_model_file
    custom.custom_output_info = json.dumps(data1)
    custom.save()


def makerThread(genome_assembly, output_dir, augustus_species, sister_proteome, num_cores,
                repeat_model_file, trinity_asms, no_genemark_fungus, custom_id):
    custom = models.Custom.objects.get(id=custom_id)
    tasks.set_loggings(os.path.join(output_dir, 'output'))
    tasks.run_maker(genome_assembly, output_dir, augustus_species, sister_proteome, num_cores,
                    repeat_model_file, trinity_asms, no_genemark_fungus)
    maker_gff3s = glob(
        os.path.join(output_dir, 'maker_out', '*/maker_*.gff3')
    )
    maker_faas = glob(os.path.join(output_dir, 'maker_out', '*/maker_*.faa'))
    data1 = {}
    data1['running_status'] = 'finished'
    data1['maker_gff3s'] = " ".join(maker_gff3s)
    data1['maker_faas'] = " ".join(maker_faas)
    custom.custom_output_info = json.dumps(data1)
    custom.save()


def repeatmaskerThread(genome_assembly, output_dir, repeat_model_file, num_cores, custom_id):
    custom = models.Custom.objects.get(id=custom_id)
    tasks.set_loggings(os.path.join(output_dir, 'output'))
    tasks.run_repeat_masker(genome_assembly, output_dir, repeat_model_file, num_cores)
    masked_assembly = glob(os.path.join(output_dir, '*.fasta.masked'))[0]
    data1 = {}
    data1['running_status'] = "finished"
    data1['masked_assembly'] = masked_assembly
    custom.custom_output_info = json.dumps(data1)
    custom.save()


def brakerThread(masked_assembly, trans_bams, output_dir, num_cores, no_braker_fungus, custom_id):
    custom = models.Custom.objects.get(id=custom_id)
    tasks.set_loggings(os.path.join(output_dir, 'output'))
    tasks.run_braker1(masked_assembly, trans_bams, output_dir, num_cores, no_braker_fungus)
    prefixes = [os.path.basename(os.path.splitext(x)[0]) for x in trans_bams]
    prefixes_u = list(set(prefixes))
    braker1_gff3s = []
    braker1_faas = []
    for prefix in prefixes_u:
        braker1_gff3 = os.path.join(
            output_dir, 'braker1_out', prefix, 'braker1_{}.gff3'.format(prefix)
        )
        braker1_gff3s.append(braker1_gff3)
        braker1_faa = os.path.join(
            output_dir, 'braker1_out', prefix, 'braker1_{}.faa'.format(prefix)
        )
        braker1_faas.append(braker1_faa)
    data1 = {}
    data1['running_status'] = 'finished'
    data1['braker1_gff3s'] = " ".join(braker1_gff3s)
    data1['braker1_faas'] = " ".join(braker1_faas)
    custom.custom_output_info = json.dumps(data1)
    custom.save()


def augustusThread(masked_assembly, output_dir, augustus_species, custom_id):
    custom = models.Custom.objects.get(id=custom_id)
    tasks.set_loggings(os.path.join(output_dir, 'output'))
    tasks.run_augustus(masked_assembly, output_dir, augustus_species)
    augustus_gff3 = os.path.join(output_dir, 'augustus.gff3')
    augustus_faa = os.path.join(output_dir, 'augustus.faa')
    data1 = {}
    data1['running_status'] = 'finished'
    data1['augustus_faa'] = augustus_faa
    data1['augustus_gff'] = augustus_gff3
    custom.custom_output_info = json.dumps(data1)
    custom.save()


def buscoThread(faa_files, output_dir, num_cores, custom_id):
    custom = models.Custom.objects.get(id=custom_id)
    tasks.set_loggings(os.path.join(output_dir, 'output'))
    tasks.run_busco_tatall(faa_files, output_dir, num_cores)
    busco_out_dir = os.path.join(output_dir, 'busco_out')
    data1 = {}
    data1['running_status'] = 'finished'
    data1['busco_out'] = busco_out_dir
    custom.custom_output_info = json.dumps(data1)
    custom.save()


def nrThread(faa_files, output_dir, custom_id):
    custom = models.Custom.objects.get(id=custom_id)
    tasks.set_loggings(os.path.join(output_dir, 'output'))
    tasks.make_nr_prot(faa_files, output_dir)
    gene_filtering_dir = os.path.join(output_dir, 'gene_filtering')
    nr_prot_file = os.path.join(gene_filtering_dir, 'nr_prot.faa')
    nr_prot_mapping_file = os.path.join(
        gene_filtering_dir, 'nr_prot_mapping.txt'
    )
    data1 = {}
    data1['running_status'] = 'finished'
    data1['nr_prot_file'] = nr_prot_file
    data1['nr_prot_mapping_file'] = nr_prot_mapping_file
    custom.custom_output_info = json.dumps(data1)
    custom.save()


def blastpThread(nr_prot_file, output_dir, sister_proteome, num_cores, custom_id):
    custom = models.Custom.objects.get(id=custom_id)
    tasks.set_loggings(os.path.join(output_dir, 'output'))
    tasks.run_blastp(nr_prot_file, output_dir, sister_proteome, num_cores)
    blastp_output = os.path.join(output_dir, 'gene_filtering', 'nr_prot.blastp')
    data1 = {}
    data1['running_status'] = 'finished'
    data1['blastp_output'] = blastp_output
    custom.custom_output_info = json.dumps(data1)
    custom.save()


def pfamscanThread(nr_prot_file, output_dir, num_cores, custom_id):
    custom = models.Custom.objects.get(id=custom_id)
    tasks.set_loggings(os.path.join(output_dir, 'output'))
    tasks.run_pfam_scan(nr_prot_file, output_dir, num_cores)
    pfam_scan_out = os.path.join(
        output_dir, 'gene_filtering', 'nr_prot.pfam_scan'
    )
    data1 = {}
    data1['running_status'] = 'finished'
    data1['pfam_scan_out'] = pfam_scan_out
    custom.custom_output_info = json.dumps(data1)
    custom.save()


def blastnThread(gff3_files, genome_assembly, trinity_asm, output_dir, custom_id):
    custom = models.Custom.objects.get(id=custom_id)
    tasks.set_loggings(os.path.join(output_dir, 'output'))
    blastn_out_files = tasks.run_blastn_totall(gff3_files, genome_assembly, trinity_asm, output_dir)
    data1 = {}
    data1['running_status'] = 'finished'
    data1['blastn_out_files'] = ' '.join(blastn_out_files)
    custom.custom_output_info = json.dumps(data1)
    custom.save()


def importblastpThread(blastp_output, nr_prot_mapping, custom_id):
    custom = models.Custom.objects.get(id=custom_id)
    tasks.set_loggings(os.path.join(blastp_output, 'output'))
    tasks.import_blastp(blastp_output, nr_prot_mapping)
    blastp_out_dir = os.path.dirname(blastp_output)
    blastp_dict = os.path.join(blastp_out_dir, 'blastp_score.p')
    data1 = {}
    data1['running_status'] = 'finished'
    data1['blastp_dict'] = blastp_dict
    custom.custom_output_info = json.dumps(data1)
    custom.save()


def importbuscoThread(busco_out_dir, custom_id):
    custom = models.Custom.objects.get(id=custom_id)
    tasks.set_loggings(os.path.join(busco_out_dir, 'output'))
    tasks.import_busco(busco_out_dir, busco_out_dir)
    gene_filtering_dir = os.path.join(busco_out_dir, 'gene_filtering')
    busco_dict = os.path.join(gene_filtering_dir, 'busco_score.p')
    data1 = {}
    data1['running_status'] = 'finished'
    data1['busco_dict'] = busco_dict
    custom.custom_output_info = json.dumps(data1)
    custom.save()


def importpfamThread(pfam_scan_out, nr_prot_mapping_file, custom_id):
    custom = models.Custom.objects.get(id=custom_id)
    tasks.set_loggings(os.path.join(pfam_scan_out, 'output'))
    tasks.import_pfam(pfam_scan_out, nr_prot_mapping_file)
    pfam_scan_out_dir = os.path.dirname(pfam_scan_out)
    pfam_dict = os.path.join(pfam_scan_out_dir, 'pfam_score.p')
    data1 = {}
    data1['running_status'] = 'finished'
    data1['pfam_dict'] = pfam_dict
    custom.custom_output_info = json.dumps(data1)
    custom.save()


def importblastnThread(blastn_output_files, custom_id):
    custom = models.Custom.objects.get(id=custom_id)
    tasks.set_loggings(os.path.join(os.path.dirname(blastn_output_files[0]), 'output'))
    output_dir = os.path.dirname(blastn_output_files[0])
    tasks.import_blastn(blastn_output_files, output_dir)
    gene_filtering_dir = os.path.join(output_dir, 'gene_filtering')
    blastn_dict = os.path.join(gene_filtering_dir, 'blastn_score.p')
    data1 = {}
    data1['running_status'] = 'finished'
    data1['blastn_dict'] = blastn_dict
    custom.custom_output_info = json.dumps(data1)
    custom.save()


def baddictThread(gff3_files, genome_assembly, output_dir, custom_id):
    custom = models.Custom.objects.get(id=custom_id)
    tasks.set_loggings(os.path.join(output_dir, 'output'))
    tasks.catch_bad_genes(gff3_files, genome_assembly, output_dir)
    bad_dict = os.path.join(os.path.join(output_dir, 'gene_filtering'), 'D_bad.p')
    data1 = {}
    data1['running_status'] = 'finished'
    data1['bad_dict'] = bad_dict
    custom.custom_output_info = json.dumps(data1)
    custom.save()


def filterThread(
        genome_assembly, gff3_files, blastp_dict, busco_dict, pfam_dict, blastn_dict,
        bad_dict, nr_prot_file, nr_prot_mapping_file, output_dir, custom_id
):
    custom = models.Custom.objects.get(id=custom_id)
    tasks.set_loggings(os.path.join(output_dir, 'output'))
    tasks.filter_gff3s(genome_assembly, gff3_files, blastp_dict, busco_dict, pfam_dict, blastn_dict,
                       bad_dict, nr_prot_file, nr_prot_mapping_file, output_dir)
    data1 = {}
    data1['running_status'] = 'finished'
    custom.custom_output_info = json.dumps(data1)
    custom.save()


def postprocessThread(genome_assembly, ori_gff3, custom_id):
    custom = models.Custom.objects.get(id=custom_id)
    tasks.set_loggings(os.path.join(ori_gff3, 'output'))
    tasks.gff3_postprocess(genome_assembly, ori_gff3)
    data1 = {}
    data1['running_status'] = 'finished'
    custom.custom_output_info = json.dumps(data1)
    custom.save()


def createmarkdownThread(genome_assembly, final_gff3, trinity_asms, trans_bams, output_dir, custom_id):
    custom = models.Custom.objects.get(id=custom_id)
    tasks.set_loggings(os.path.join(output_dir, 'output'))
    tasks.create_markdown(genome_assembly, output_dir, trans_bams, trinity_asms, final_gff3)
    data1 = {}
    data1['running_status'] = 'finished'
    custom.custom_output_info = json.dumps(data1)
    custom.save()


def delete_custom(request):
    custom_id = request.GET.get('custom_id')
    custom = models.Custom.objects.get(id=custom_id)
    running_info = json.loads(custom.custom_output_info)
    if running_info['running_status'] == 'finished':
        custom.delete()
        return HttpResponse('success')
    else:
        return HttpResponse('')


def delete_project(request):
    project_name = request.GET.get('project_name')
    project = models.Project.objects.get(project_name=project_name)
    status = project.status
    if status == 'Running':
        return HttpResponse('')
    else:
        project.delete()
        return HttpResponse('success')


def get_markdown(request):
    project_name = request.GET.get('project_name')
    markdown_path = os.path.join(BASE_DIR, 'Ace', project_name, 'output', 'fungap_out', 'fungap_out.html')
    with open(markdown_path, 'r') as f:
        content = f.read()
    pattern = re.compile(r'src=".*"')
    figure1_path = '/static/fungap/images/' + project_name + '/fungap_out_trans_len_dist.png'
    figure2_path = '/static/fungap/images/' + project_name + '/fungap_out_prot_len_dist.png'
    content = pattern.sub("src='{}' style='width: 100%; height: auto'".format(figure1_path), content, 1)
    content = pattern.sub('src="{}" style="width: 100%; height: auto"'.format(figure2_path), content, 1)
    patter1 = re.compile(r'<head>[\s\S]*</head>')
    content = patter1.sub('', content, 1)
    content = content.replace('<body>', '')
    content = content.replace('</body>', '')
    return HttpResponse(content)


def get_log(request):
    project_name = request.GET.get('project_name')
    log_path = os.path.join(BASE_DIR, 'Ace', project_name, 'output', 'logs', 'fungap.log')
    with open(log_path, 'r') as f:
        content = f.read()
    content = content.replace('\n', '<br>')
    return HttpResponse(content)
