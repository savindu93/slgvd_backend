from django.db import models
from django.contrib.auth.models import AbstractBaseUser, BaseUserManager, PermissionsMixin

# The CustomUserManager has custom defined create_user and create_superuser functions. The custom create_user expects the and password during user creation instead of the default username.
class CustomUserManager(BaseUserManager):

    def create_user(self, email, password, **extra_fields):

        if not email:
            raise ValueError("The Email field must be required")
        
        email = self.normalize_email(email)
        user = self.model(email = email, **extra_fields)
        user.set_password(password)
        user.save(using=self.db)
        return user
    
    def create_superuser(self, email, password=None, **extra_fields):

        extra_fields.setdefault('is_staff', True)
        extra_fields.setdefault('is_superuser', True)

        return self.create_user(email, password, **extra_fields)
    
class CustomUser(AbstractBaseUser, PermissionsMixin):

    username = models.CharField(max_length = 150, unique=True)
    email = models.EmailField(unique=True)
    first_name = models.CharField(max_length=150)
    last_name = models.CharField(max_length=150)
    is_staff = models.BooleanField(default=False)
    is_active = models.BooleanField(default=True)
    is_superuser = models.BooleanField(default=False)
    date_joined = models.DateTimeField(auto_now_add=True)

    USERNAME_FIELD = 'email'
    REQUIRED_FIELDS = []

    objects = CustomUserManager()

class Submission(models.Model):

    submission_id = models.AutoField(unique = True, primary_key=True)
    no_individuals = models.IntegerField()
    submission_date = models.DateField(auto_now_add=True)
    username = models.ForeignKey(CustomUser, to_field = 'username',on_delete=models.CASCADE)

    class Meta:
        db_table = 'submissions'

class Variant(models.Model):

    variation_id = models.CharField(max_length = 255, primary_key = True)
    chromosome = models.CharField(max_length = 8)
    position = models.IntegerField()
    ref_allele = models.CharField(max_length = 255)
    alt_allele = models.CharField(max_length = 255)
    gene_name = models.CharField(max_length = 255)
    consequence = models.CharField(max_length=255, null=True)
    submission = models.ForeignKey(Submission, on_delete=models.CASCADE)

    class Meta:
         db_table = 'ls50_variants'

class CNV(models.Model):

    variation_id = models.CharField(max_length = 255, primary_key = True)
    chromosome = models.CharField(max_length = 8)
    start_pos = models.IntegerField()
    end_pos = models.IntegerField()
    site_count = models.IntegerField()
    consequence = models.CharField(max_length=255, null=True)
    submission = models.ForeignKey(Submission, on_delete=models.CASCADE)
    last_updated = models.JSONField(null = True)


    class Meta:
        db_table = 'mt50_variants'

class Frequency(models.Model):

    variation = models.OneToOneField(Variant, to_field='variation_id', on_delete = models.CASCADE, primary_key=True)
    homo_count = models.IntegerField(default=0)
    het_count = models.IntegerField(default=0)
    last_updated = models.JSONField(null = True)

    class Meta:

        db_table = 'allele_count'

# Takes the total count of variants by type of variant and type of cosequence of each variant 
class VarCounts(models.Model):

    entry_id = models.AutoField(unique = True, primary_key = True)
    cnv_count = models.IntegerField(default = 0)
    ssv_count = models.IntegerField(default = 0)
    exo_count = models.IntegerField(default = 0)
    intro_count = models.IntegerField(default = 0)
    inter_count = models.IntegerField(default = 0)
    reg_count = models.IntegerField(default = 0)
    spli_count = models.IntegerField(default = 0)
    intra_count = models.IntegerField(default = 0)

    class Meta:
        db_table = 'var_counts'






